#include "params.h"
#include <cuda_runtime.h>
#include <iostream>
#include <chrono>
#include <vector>

#if N < 65536u
typedef uint16_t MannaSizeType;
#else
typedef uint32_t MannaSizeType;
#endif

typedef uint32_t MannaItemType;

__device__ inline uint32_t xorshift32(uint32_t &state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

// GPU-based reduction to check if any block is active
__global__ void check_activity_kernel(uint32_t *active_mask, int num_blocks, int *global_active) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Use shared memory for block-level reduction
    __shared__ int block_has_active[256]; // Assuming max 256 threads per block
    
    // Each thread checks its assigned blocks
    bool thread_active = false;
    for (int i = idx; i < num_blocks; i += blockDim.x * gridDim.x) {
        if (active_mask[i] != 0) {
            thread_active = true;
            break;
        }
    }
    
    block_has_active[threadIdx.x] = thread_active ? 1 : 0;
    __syncthreads();
    
    // Block-level reduction
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            block_has_active[threadIdx.x] = block_has_active[threadIdx.x] || block_has_active[threadIdx.x + stride];
        }
        __syncthreads();
    }
    
    // First thread of each block writes to global memory
    if (threadIdx.x == 0 && block_has_active[0]) {
        atomicOr(global_active, 1);
    }
}

// Optimized destabilization kernel with better memory access patterns
__global__ void destab_kernel_opt(
    MannaItemType *h,
    MannaItemType *temp_h,
    uint32_t      *rng)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= N) return;

    MannaItemType c = h[gid];
    
    // Early exit for inactive cells
    if (c <= 1) return;
    
    // Handle overflow case
    if (c >= 32) {
        h[gid] = 0;
        return;
    }

    h[gid] = 0;
    uint32_t r = xorshift32(rng[gid]);
    uint32_t mask = ((1u << c) - 1u) & r;
    uint32_t to_right = __popc(mask);
    uint32_t to_left  = c - to_right;

    unsigned right = (gid + 1) & (N - 1);
    unsigned left  = (gid + N - 1) & (N - 1);

    // Use non-atomic writes when possible (if we can guarantee no conflicts)
    atomicAdd(&temp_h[right], to_right);
    atomicAdd(&temp_h[left],  to_left);
}

// Optimized merge kernel with early termination
__global__ void merge_kernel_opt(
    MannaItemType *h,
    MannaItemType *temp_h,
    uint32_t      *active_mask)
{
    unsigned int gid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int lane = threadIdx.x & 31;

    MannaItemType new_val = 0;
    bool my_active = false;
    
    // Merge step
    if (gid < N) {
        new_val = h[gid] + temp_h[gid];
        h[gid] = new_val;
        my_active = (new_val > 1);
        temp_h[gid] = 0; // Clear for next iteration
    }

    // Warp-wide ballot for activity detection
    uint32_t warp_mask = __ballot_sync(0xFFFFFFFFu, my_active);

    // Use shared memory to avoid redundant atomic operations
    __shared__ bool block_active;
    if (threadIdx.x == 0) {
        block_active = false;
    }
    __syncthreads();
    
    if (lane == 0 && warp_mask != 0) {
        block_active = true;
    }
    __syncthreads();

    // Write block activity once per block
    if (threadIdx.x == 0) {
        active_mask[blockIdx.x] = block_active ? 1u : 0u;
    }
}

// Batch multiple iterations in a single kernel launch
__global__ void multi_step_kernel(
    MannaItemType *h,
    MannaItemType *temp_h,
    uint32_t      *rng_state,
    uint32_t      *active_mask,
    int           max_steps,
    int           *steps_taken)
{
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    
    for (int step = 0; step < max_steps; step++) {
        __syncthreads(); // Synchronize across all threads in grid (requires cooperative groups for multi-block)
        
        // Check if we should continue (simplified version)
        if (gid == 0 && step > 0) {
            // Simple activity check - in practice you'd need a more sophisticated approach
            bool any_active = false;
            for (int i = 0; i < gridDim.x && !any_active; i++) {
                if (active_mask[i]) any_active = true;
            }
            if (!any_active) {
                *steps_taken = step;
                return;
            }
        }
        
        // Reset temp array
        if (gid < N) {
            temp_h[gid] = 0;
        }
        if (gid < gridDim.x) {
            active_mask[gid] = 0;
        }
        __syncthreads();
        
        // Destabilization phase
        if (gid < N) {
            MannaItemType c = h[gid];
            if (c > 1 && c < 32) {
                h[gid] = 0;
                uint32_t r = xorshift32(rng_state[gid]);
                uint32_t mask = ((1u << c) - 1u) & r;
                uint32_t to_right = __popc(mask);
                uint32_t to_left = c - to_right;
                
                unsigned right = (gid + 1) & (N - 1);
                unsigned left = (gid + N - 1) & (N - 1);
                
                atomicAdd(&temp_h[right], to_right);
                atomicAdd(&temp_h[left], to_left);
            } else if (c >= 32) {
                h[gid] = 0;
            }
        }
        __syncthreads();
        
        // Merge phase
        bool my_active = false;
        if (gid < N) {
            h[gid] += temp_h[gid];
            my_active = (h[gid] > 1);
        }
        
        // Activity detection
        unsigned int lane = threadIdx.x & 31;
        uint32_t warp_mask = __ballot_sync(0xFFFFFFFFu, my_active);
        
        __shared__ bool block_active;
        if (threadIdx.x == 0) {
            block_active = false;
        }
        __syncthreads();
        
        if (lane == 0 && warp_mask != 0) {
            block_active = true;
        }
        __syncthreads();
        
        if (threadIdx.x == 0) {
            active_mask[blockIdx.x] = block_active ? 1u : 0u;
        }
    }
    
    if (gid == 0) {
        *steps_taken = max_steps;
    }
}

// Other kernels remain the same
__global__ void inicializacion_kernel(MannaItemType *h) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        h[idx] = static_cast<MannaItemType>((idx + 1) * DENSITY)
               - static_cast<MannaItemType>(idx * DENSITY);
    }
}

__global__ void desestabilizacion_kernel(
    MannaItemType *h, MannaItemType *temp_h, uint32_t *rng)
{
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N) return;
    if (h[idx] == 1) {
        bool dir = xorshift32(rng[idx]) & 1;
        int j = (idx + 2*dir - 1) & (N-1);
        atomicAdd(&temp_h[j], 1u);
        h[idx] = 0;
    }
}

__global__ void setup_rng(uint32_t *rng_states, uint32_t seed) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        rng_states[idx] = seed ^ (idx * 0x9E3779B1u);
    }
}

class MannaCUDA {
private:
    MannaItemType *h, *temp_h;
    uint32_t      *rng_state;
    uint32_t      *active_mask;
    int           *global_active; // GPU-based activity flag
    int           *steps_counter; // GPU-based step counter
    int            block_size, grid_size;

public:
    MannaCUDA(){
        block_size = BLOCK_SIZE < N ? BLOCK_SIZE : N;
        grid_size  = (N + block_size - 1)/ block_size;

        cudaMalloc(&h,            N * sizeof(*h));
        cudaMalloc(&temp_h,       N * sizeof(*temp_h));
        cudaMalloc(&rng_state,    N * sizeof(*rng_state));
        cudaMalloc(&active_mask,  grid_size * sizeof(*active_mask));
        cudaMalloc(&global_active, sizeof(int));
        cudaMalloc(&steps_counter, sizeof(int));

        setup_rng<<<grid_size,block_size>>>(rng_state, SEED);
        cudaDeviceSynchronize();
    }
    
    ~MannaCUDA(){
        cudaFree(h);
        cudaFree(temp_h);
        cudaFree(rng_state);
        cudaFree(active_mask);
        cudaFree(global_active);
        cudaFree(steps_counter);
    }

    void inicializacion(){
        inicializacion_kernel<<<grid_size, block_size>>>(h);
        cudaDeviceSynchronize();
    }

    void desestabilizacion_inicial(){
        cudaMemsetAsync(temp_h, 0, N * sizeof(*temp_h));
        desestabilizacion_kernel<<<grid_size, block_size>>>(h, temp_h, rng_state);
        cudaDeviceSynchronize();
    }

    // Optimized version that avoids host-device transfers
    bool descargar_opt() {
        // Clear arrays (do this less frequently or batch)
        cudaMemsetAsync(temp_h, 0, N * sizeof(*temp_h));
        cudaMemsetAsync(active_mask, 0, grid_size * sizeof(*active_mask));
        cudaMemsetAsync(global_active, 0, sizeof(int));

        // Launch kernels
        destab_kernel_opt<<<grid_size, block_size>>>(h, temp_h, rng_state);
        merge_kernel_opt<<<grid_size, block_size>>>(h, temp_h, active_mask);

        // Check activity on GPU
        int check_blocks = (grid_size + 255) / 256;
        check_activity_kernel<<<check_blocks, 256>>>(active_mask, grid_size, global_active);

        // Only one small memory transfer
        int host_active;
        cudaMemcpy(&host_active, global_active, sizeof(int), cudaMemcpyDeviceToHost);
        
        return host_active != 0;
    }

    // Original version for comparison
    bool descargar() {
        cudaMemsetAsync(temp_h, 0, N * sizeof(*temp_h));
        cudaMemsetAsync(active_mask, 0, grid_size * sizeof(*active_mask));

        destab_kernel_opt<<<grid_size, block_size>>>(h, temp_h, rng_state);
        merge_kernel_opt<<<grid_size, block_size>>>(h, temp_h, active_mask);

        std::vector<uint32_t> host_mask(grid_size);
        cudaMemcpy(host_mask.data(), active_mask, 
                   grid_size * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        
        for (auto &w : host_mask) {
            if (w) return true;
        }
        return false;
    }

    // Batch processing version
    int run_batched_steps(int max_batch_size) {
        cudaMemset(steps_counter, 0, sizeof(int));
        
        // This is a simplified version - full implementation would need cooperative groups
        // for proper multi-block synchronization
        multi_step_kernel<<<grid_size, block_size>>>(
            h, temp_h, rng_state, active_mask, max_batch_size, steps_counter);
        
        int steps_taken;
        cudaMemcpy(&steps_taken, steps_counter, sizeof(int), cudaMemcpyDeviceToHost);
        return steps_taken;
    }

    void swap_arrays() {
        MannaItemType *temp = h;
        h = temp_h;
        temp_h = temp;
    }

    void print_array_h() {
        std::vector<MannaItemType> host_h(N);
        cudaMemcpy(host_h.data(), h, N * sizeof(MannaItemType), cudaMemcpyDeviceToHost);
        for (size_t i = 0; i < N; ++i) {
            std::cout << host_h[i] << " ";
        }
        std::cout << "\n";
    }
};

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    MannaCUDA manna;

    manna.inicializacion();
    manna.desestabilizacion_inicial();
    manna.swap_arrays();

    uint32_t t = 0;
    bool active;
    
    do {
        active = manna.descargar_opt(); // Use optimized version
        ++t;
    } while (active && t < NSTEPS);

    uint32_t processed = t * N * DENSITY;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "=== RESULTADOS FINALES ===\n";
    std::cout << "Steps taken: " << t << "\n";
    std::cout << "Tiempo (s): " << static_cast<double>(duration.count())/1e6 << "\n";
    std::cout << "Granos procesados: " << processed << "\n";
    std::cout << "Granos/us: " << static_cast<double>(processed)/duration.count() << "\n";
    
    return 0;
}
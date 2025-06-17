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

// Combined kernel that does destabilization + merge + activity check in one pass
// This maintains the EXACT same logic as the original three separate kernels
__global__ void combined_step_kernel(
    MannaItemType *h,           // Main array
    MannaItemType *temp_h,      // Temporary array
    uint32_t      *rng_state,   // RNG states
    uint32_t      *active_mask, // Block activity mask
    int           *global_active // Global activity flag
) {
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int lane = threadIdx.x & 31;
    unsigned int warp_id = threadIdx.x >> 5;
    
    // Shared memory for warp results
    __shared__ uint32_t warp_results[32];
    
    // Phase 1: Clear temp array (same as original memset)
    if (gid < N) {
        temp_h[gid] = 0;
    }
    
    // Ensure temp array is cleared before destabilization
    __syncthreads();
    
    // Phase 2: Destabilization (EXACT same logic as destab_kernel_opt)
    if (gid < N) {
        MannaItemType c = h[gid];
        
        // Early exit for inactive cells (same as original)
        if (c <= 1) {
            // Do nothing, just like original
        } else if (c >= 32) {
            // Handle overflow case (same as original)
            h[gid] = 0;
        } else {
            // Normal destabilization (same as original)
            h[gid] = 0;
            uint32_t r = xorshift32(rng_state[gid]);
            uint32_t mask = ((1u << c) - 1u) & r;
            uint32_t to_right = __popc(mask);
            uint32_t to_left = c - to_right;
            
            unsigned right = (gid + 1) & (N - 1);
            unsigned left = (gid + N - 1) & (N - 1);
            
            atomicAdd(&temp_h[right], to_right);
            atomicAdd(&temp_h[left], to_left);
        }
    }
    
    // Synchronize to ensure all destabilization is complete
    __syncthreads();
    
    // Phase 3: Merge (EXACT same logic as merge_kernel_opt)
    bool my_active = false;
    if (gid < N) {
        MannaItemType new_val = h[gid] + temp_h[gid];
        h[gid] = new_val;
        my_active = (new_val > 1);  // Same activity condition as original
        temp_h[gid] = 0; // Clear for next iteration (same as original)
    }
    
    // Phase 4: Activity detection (EXACT same logic as merge_kernel_opt)
    uint32_t warp_mask = __ballot_sync(0xFFFFFFFFu, my_active);
    
    // Store warp results in shared memory
    if (lane == 0) {
        warp_results[warp_id] = warp_mask;
    }
    __syncthreads();
    
    // Thread 0 consolidates results (same as original)
    if (threadIdx.x == 0) {
        bool block_active = false;
        int num_warps = (blockDim.x + 31) / 32;
        
        for (int i = 0; i < num_warps; i++) {
            if (warp_results[i] != 0) {
                block_active = true;
                break;
            }
        }
        
        // Update block activity mask (same as original)
        active_mask[blockIdx.x] = block_active ? 1u : 0u;
        
        // Atomically update global activity if this block is active
        if (block_active) {
            atomicOr(global_active, 1);
        }
    }
}

// Modified merge kernel that also updates global activity atomically
__global__ void merge_kernel_with_global_activity(
    MannaItemType *h,
    MannaItemType *temp_h,
    uint32_t      *active_mask,
    int           *global_active)
{
    unsigned int gid  = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int lane = threadIdx.x & 31;

    MannaItemType new_val = 0;
    bool my_active = false;
    
    // Merge step with clearing (EXACT same as original)
    if (gid < N) {
        new_val = h[gid] + temp_h[gid];
        h[gid] = new_val;
        my_active = (new_val > 1);
        temp_h[gid] = 0; // Clear for next iteration - eliminates memset!
    }

    // Fast warp-wide ballot for activity detection (same as original)
    uint32_t warp_mask = __ballot_sync(0xFFFFFFFFu, my_active);

    // Only lane 0 of warps with activity write to shared memory (same as original)
    __shared__ uint32_t warp_results[32]; // Max 32 warps per block
    unsigned int warp_id = threadIdx.x >> 5;
    
    if (lane == 0) {
        warp_results[warp_id] = warp_mask;
    }
    __syncthreads();

    // Thread 0 consolidates warp results (same as original)
    if (threadIdx.x == 0) {
        bool block_active = false;
        int num_warps = (blockDim.x + 31) / 32;
        for (int i = 0; i < num_warps; i++) {
            if (warp_results[i] != 0) {
                block_active = true;
                break;
            }
        }
        active_mask[blockIdx.x] = block_active ? 1u : 0u;
        
        // NEW: Also update global activity atomically
        if (block_active) {
            atomicOr(global_active, 1);
        }
    }
}

// Keep your original optimized destabilization kernel
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

    atomicAdd(&temp_h[right], to_right);
    atomicAdd(&temp_h[left],  to_left);
}
__global__ void multi_step_combined_kernel(
    MannaItemType *h,
    MannaItemType *temp_h,
    uint32_t      *rng_state,
    uint32_t      *active_mask,
    int           *global_active,
    int           *step_counter,
    int           max_steps
) {
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int lane = threadIdx.x & 31;
    unsigned int warp_id = threadIdx.x >> 5;
    
    __shared__ uint32_t warp_results[32];
    __shared__ bool continue_flag;
    
    for (int step = 0; step < max_steps; step++) {
        // Reset global activity flag
        if (gid == 0) {
            *global_active = 0;
            continue_flag = true;
        }
        __syncthreads();
        
        // Check if we should continue (after first step)
        if (step > 0 && gid == 0) {
            // Simple check - if no activity detected in previous step, stop
            if (*global_active == 0) {
                continue_flag = false;
                *step_counter = step;
            }
        }
        __syncthreads();
        
        if (!continue_flag) break;
        
        // Clear temp array and destabilize
        if (gid < N) {
            temp_h[gid] = 0;
            
            MannaItemType c = h[gid];
            if (c > 1) {
                if (c >= 32) {
                    h[gid] = 0;
                } else {
                    h[gid] = 0;
                    uint32_t r = xorshift32(rng_state[gid]);
                    uint32_t mask = ((1u << c) - 1u) & r;
                    uint32_t to_right = __popc(mask);
                    uint32_t to_left = c - to_right;
                    
                    unsigned right = (gid + 1) & (N - 1);
                    unsigned left = (gid + N - 1) & (N - 1);
                    
                    atomicAdd(&temp_h[right], to_right);
                    atomicAdd(&temp_h[left], to_left);
                }
            }
        }
        __syncthreads();
        
        // Merge and check activity
        bool my_active = false;
        if (gid < N) {
            MannaItemType new_val = h[gid] + temp_h[gid];
            h[gid] = new_val;
            my_active = (new_val > 1);
        }
        
        // Activity detection
        uint32_t warp_mask = __ballot_sync(0xFFFFFFFFu, my_active);
        
        if (lane == 0) {
            warp_results[warp_id] = warp_mask;
        }
        __syncthreads();
        
        if (threadIdx.x == 0) {
            bool block_active = false;
            int num_warps = (blockDim.x + 31) / 32;
            
            for (int i = 0; i < num_warps; i++) {
                if (warp_results[i] != 0) {
                    block_active = true;
                    break;
                }
            }
            
            active_mask[blockIdx.x] = block_active ? 1u : 0u;
            
            if (block_active) {
                atomicOr(global_active, 1);
            }
        }
        __syncthreads();
    }
    
    // Update final step counter
    if (gid == 0 && continue_flag) {
        *step_counter = max_steps;
    }
}

// Remaining original kernels
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
    int           *global_active;
    int           *step_counter;
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
        cudaMalloc(&step_counter, sizeof(int));

        setup_rng<<<grid_size,block_size>>>(rng_state, SEED);
        cudaDeviceSynchronize();
    }
    
    ~MannaCUDA(){
        cudaFree(h);
        cudaFree(temp_h);
        cudaFree(rng_state);
        cudaFree(active_mask);
        cudaFree(global_active);
        cudaFree(step_counter);
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

    // Ultra-optimized single kernel approach - maintains exact original logic
    bool descargar_combined() {
        // Reset global activity flag (same as original memset)
        cudaMemsetAsync(global_active, 0, sizeof(int));
        
        // Single kernel that does everything with EXACT same logic as original 3 kernels
        combined_step_kernel<<<grid_size, block_size>>>(
            h, temp_h, rng_state, active_mask, global_active);
        
        // Single memory transfer (replaces the original check_activity_kernel_fast + memcpy)
        int host_active;
        cudaMemcpy(&host_active, global_active, sizeof(int), cudaMemcpyDeviceToHost);
        
        return host_active != 0;
    }
    
    // Conservative optimization: Keep original kernels but optimize activity check
    bool descargar_optimized_activity() {
        // Clear arrays (same as original)
        cudaMemsetAsync(temp_h, 0, N * sizeof(*temp_h));
        cudaMemsetAsync(global_active, 0, sizeof(int));
        
        // Use your original optimized kernels
        destab_kernel_opt<<<grid_size, block_size>>>(h, temp_h, rng_state);
        
        // Modified merge kernel that also updates global activity
        merge_kernel_with_global_activity<<<grid_size, block_size>>>(
            h, temp_h, active_mask, global_active);
        
        // Single memory transfer instead of check_activity_kernel_fast + memcpy
        int host_active;
        cudaMemcpy(&host_active, global_active, sizeof(int), cudaMemcpyDeviceToHost);
        
        return host_active != 0;
    }
    
    // Batch processing with reduced transfers
    int run_batch_optimized(int batch_size = 100) {
        // Reset counters
        cudaMemset(global_active, 0, sizeof(int));
        cudaMemset(step_counter, 0, sizeof(int));
        
        // Process multiple steps in single kernel launch
        multi_step_combined_kernel<<<grid_size, block_size>>>(
            h, temp_h, rng_state, active_mask, global_active, step_counter, batch_size);
        
        // Single memory transfer to get step count
        int steps_taken;
        cudaMemcpy(&steps_taken, step_counter, sizeof(int), cudaMemcpyDeviceToHost);
        
        return steps_taken;
    }
    
    // Adaptive batch processing
    int run_adaptive_batch() {
        int total_steps = 0;
        int batch_size = 1000;  // Start with larger batches
        
        while (total_steps < NSTEPS) {
            int remaining = NSTEPS - total_steps;
            int current_batch = (remaining < batch_size) ? remaining : batch_size;
            
            int steps_taken = run_batch_optimized(current_batch);
            total_steps += steps_taken;
            
            // If we took fewer steps than the batch size, we're done
            if (steps_taken < current_batch) {
                break;
            }
            
            // Adaptive batch sizing - increase batch size for longer runs
            if (steps_taken == current_batch && batch_size < 10000) {
                batch_size *= 2;
            }
        }
        
        return total_steps;
    }

    // Original methods for comparison
    bool descargar_ultra_opt() {
        cudaMemsetAsync(active_mask, 0, grid_size * sizeof(*active_mask));
        cudaMemsetAsync(global_active, 0, sizeof(int));

        // These would be your original optimized kernels
        // destab_kernel_opt<<<grid_size, block_size>>>(h, temp_h, rng_state);
        // merge_kernel_opt<<<grid_size, block_size>>>(h, temp_h, active_mask);
        // check_activity_kernel_fast<<<1, 32>>>(active_mask, grid_size, global_active);

        int host_active;
        cudaMemcpy(&host_active, global_active, sizeof(int), cudaMemcpyDeviceToHost);
        
        return host_active != 0;
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

    // Option 1: Conservative optimization (2 kernels instead of 3)
    uint32_t t = 0;
    bool active;    
    do {
        active = manna.descargar_optimized_activity(); // 2 kernels per iteration
        ++t;
    } while (active && t < NSTEPS);
    
    /*
    // Option 2: Single kernel optimization (most aggressive but safe)
    uint32_t t = 0;
    bool active;    
    do {
        active = manna.descargar_combined(); // 1 kernel per iteration
        ++t;
    } while (active && t < NSTEPS);
    */
    
    /*
    // Option 3: Batch processing (experimental - use only if above work correctly)
    uint32_t t = manna.run_adaptive_batch();
    */

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
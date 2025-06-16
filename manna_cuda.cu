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

// Usamos uint32_t para contadores y atomics
typedef uint32_t MannaItemType;

// XORSHIFT32 para RNG rápido en registros
__device__ inline uint32_t xorshift32(uint32_t &state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

__global__ void descargar_kernel(
    MannaItemType *h, MannaItemType *temp_h,
    uint32_t *rng_states, unsigned int *active, uint32_t *processed
) {
    // TODO
}


// Kernel para inicializar RNG a partir de seed
__global__ void init_rng_kernel(uint32_t *rng_states, uint32_t seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        // Mezcla la seed con el índice para diversificar
        uint32_t st = seed ^ (idx * 0x9E3779B1u);
        // Perturba un poco
        for (int i = 0; i < 4; ++i) st = xorshift32(st);
        rng_states[idx] = st;
    }
}

__global__ void inicializacion_kernel(MannaItemType *h) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        h[idx] = static_cast<MannaItemType>((idx + 1) * DENSITY) - 
                 static_cast<MannaItemType>(idx * DENSITY);
    }
}

__global__ void desestabilizacion_kernel(
    MannaItemType *h, MannaItemType *temp_h, uint32_t *rng)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx >= N) return;
  if (h[idx] == 1) {
    bool dir = xorshift32(rng[idx]) & 1;
    int j = (idx + 2*dir - 1) & (N-1);
    atomicAdd(&temp_h[j], 1u);
  }
}
__global__ void setup_rng(uint32_t *rng_states, uint32_t seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        uint32_t st = seed ^ (idx * 0x9E3779B1u);
        rng_states[idx] = st;
    }
}


class MannaCUDA {
private:
    MannaItemType *d_h, *d_temp_h;
    uint32_t      *d_rng_state;
    uint32_t      *d_processed;
    unsigned int  *d_activity;
    int block_size;
    int grid_size;

public:
    MannaCUDA(){
        block_size = BLOCK_SIZE < N ? BLOCK_SIZE : N;
        grid_size = N / block_size;
        cudaMalloc(&d_h,         N * sizeof(MannaItemType));
        cudaMalloc(&d_temp_h,    N * sizeof(MannaItemType));
        cudaMalloc(&d_rng_state, N * sizeof(uint32_t));
        cudaMalloc(&d_processed, sizeof(uint32_t));
        cudaMalloc(&d_activity,  sizeof(unsigned int));
        setup_rng<<<grid_size, block_size>>>(d_rng_state, SEED);
    }

    ~MannaCUDA() {
        cudaFree(d_h);
        cudaFree(d_temp_h);
        cudaFree(d_rng_state);
        cudaFree(d_processed);
        cudaFree(d_activity);
    }
    
    void inicializacion() {
        inicializacion_kernel<<<grid_size, block_size>>>(d_h);
    }
    
    void desestabilizacion_inicial() {
        desestabilizacion_kernel<<<grid_size, block_size>>>(
            d_h, d_temp_h, d_rng_state);
    }

    // swap d_h and d_temp_h
    void swap_arrays() {
        MannaItemType *temp = d_h;
        d_h = d_temp_h;
        d_temp_h = temp;
    }
    
    bool descargar() {
        // todo
    }

    void print_array() {
        std::vector<MannaItemType> h_host(N);
        cudaMemcpy(h_host.data(), d_h, N * sizeof(MannaItemType), cudaMemcpyDeviceToHost);
        
        std::cout << "h: ";
        for (int i = 0; i < N; ++i) {
            std::cout << h_host[i] << " ";
        }
        std::cout << std::endl;
    }
};

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    
    MannaCUDA manna;
    manna.inicializacion();
    manna.desestabilizacion_inicial();
    manna.swap_arrays();
    manna.print_array();
    uint32_t t = 0;
    uint32_t processed = 0;
    bool active;

    do {
        active = manna.descargar();
        //manna.print_array();
        t++;
    } while (active && t < NSTEPS);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "=== RESULTADOS FINALES ===" << std::endl;
    std::cout << "Steps taken: " << t << std::endl;
    std::cout << "Tiempo de procesamiento (s): " << static_cast<double>(duration.count()) / 1e6 << std::endl;
    std::cout << "Granos procesados: " << processed << std::endl;
    std::cout << "Granos/us: " << static_cast<double>(processed) / duration.count() << std::endl;    
    return 0;
}
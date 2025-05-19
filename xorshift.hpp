#pragma once
#include "params.h"
#include <cstdint>
#include <immintrin.h>  // AVX2
#include <omp.h>

static constexpr int LANES = 8;
static thread_local __m256i xrng_state256;

static thread_local uint32_t bit_buffer;
static thread_local uint8_t  bit_count;


inline void xrng256_init(uint32_t seed_base, int tid) {
    const uint32_t KNUTH = 0x9e3779b9u;
    uint32_t thread_stride = KNUTH * LANES * THREADS;
    alignas(32) uint32_t seeds[LANES];
    for(int i = 0; i < LANES; ++i) {
        seeds[i] = seed_base + thread_stride * tid + KNUTH * i;
    }
    xrng_state256 = _mm256_load_si256((__m256i*)seeds);
    bit_buffer = 0;
    bit_count  = 0;
}

// Genera un vector de 256 bits aleatorios
#pragma omp declare simd
inline __m256i xrng256_next() {
    __m256i x = xrng_state256;
    __m256i t = _mm256_slli_epi32(x, 13);
    x = _mm256_xor_si256(x, t);
    t = _mm256_srli_epi32(x, 17);
    x = _mm256_xor_si256(x, t);
    t = _mm256_slli_epi32(x, 5);
    x = _mm256_xor_si256(x, t);
    xrng_state256 = x;
    return x;
}

// Extrae un unico bit aleatorio como bool
#pragma omp declare simd
inline bool xrng_next_bit() {
    if (bit_count == 0) {
        __m256i v = xrng256_next();
        uint32_t word = _mm256_extract_epi32(v, 0);
        bit_buffer = word;
        bit_count  = 32;
    }
    bool b = (bit_buffer & 1u);
    bit_buffer >>= 1;
    --bit_count;
    return b;
}
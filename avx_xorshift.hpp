// xorshift.hpp
#pragma once
#include <cstdint>
#include <immintrin.h> // AVX2

static __m256i xrng_state256;

// Inicializacion: 8 semillas de 32 bits
inline void xrng256_init(uint32_t seed_base) {
    uint32_t seeds[8];
    for(int i = 0; i < 8; ++i)
        seeds[i] = seed_base + 0x9e3779b9u * i;
    xrng_state256 = _mm256_loadu_si256((__m256i*)seeds);
}

// Bloque de 8×32 bits de RNG
inline __m256i xrng256_next() {
    __m256i x = xrng_state256;
    __m256i t;
    t = _mm256_slli_epi32(x, 13); x = _mm256_xor_si256(x, t);
    t = _mm256_srli_epi32(x, 17); x = _mm256_xor_si256(x, t);
    t = _mm256_slli_epi32(x, 5);  x = _mm256_xor_si256(x, t);
    xrng_state256 = x;
    return x;
}

static uint32_t bit_buffer[8];
static uint8_t  bit_count[8];

inline bool xrng_next_bit() {
    if (--bit_count[0] == 0) {
        alignas(32) uint32_t buf[8];
        _mm256_store_si256((__m256i*)buf, xrng256_next());
        bit_buffer[0] = buf[0];
        bit_count[0]  = 32;
    }
    bool b = bit_buffer[0] & 1u;
    bit_buffer[0] >>= 1;
    return b;
}
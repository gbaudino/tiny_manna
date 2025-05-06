// xorshift.hpp
#pragma once
#include <cstdint>
#include <immintrin.h> // AVX2
#include <cstring>


static __m256i xrng_state256;

// Inicializacion: 8 semillas de 32 bits
inline void xrng256_init(uint32_t seed_base) {
    uint32_t seeds[8];
    for(int i=0;i<8;++i) seeds[i] = seed_base + 0x9e3779b9u * i;
    xrng_state256 = _mm256_loadu_si256((__m256i*)seeds);
}

inline __m256i xrng256_next() {
    __m256i x = xrng_state256;
    // x ^= x << 13;
    __m256i t = _mm256_slli_epi32(x, 13);
    x = _mm256_xor_si256(x, t);
    // x ^= x >> 17;
    t = _mm256_srli_epi32(x, 17);
    x = _mm256_xor_si256(x, t);
    // x ^= x << 5;
    t = _mm256_slli_epi32(x, 5);
    x = _mm256_xor_si256(x, t);
    xrng_state256 = x;
    return x;
}

static uint32_t bit_buffer[8];
static uint8_t  bit_count[8];

inline bool xrng_next_bit_idx(int lane) {
    if (--bit_count[lane] == 0) {
        alignas(32) uint32_t buf[8];
        _mm256_store_si256((__m256i*)buf, xrng256_next());
        bit_buffer[lane] = buf[lane];
        bit_count[lane]  = 32;
    }
    bool b = bit_buffer[lane] & 1u;
    bit_buffer[lane] >>= 1;
    return b;
}
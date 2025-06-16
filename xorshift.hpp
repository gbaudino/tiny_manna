#pragma once
#include <cstdint>

static uint32_t xrng_state;
static uint32_t bit_buffer;
static uint8_t  bit_count;

inline void xrng_init(uint32_t seed) {
    xrng_state = seed;
    bit_buffer = 0;
    bit_count  = 0;
}

inline uint32_t xrng_next() {
    uint32_t x = xrng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    xrng_state = x;
    return x;
}

inline bool xrng_bit() {
    if (bit_count == 0) {
        bit_buffer = xrng_next();
        bit_count  = 32;
    }
    bool b = (bit_buffer & 1u);
    bit_buffer >>= 1;
    --bit_count;
    return b;
}
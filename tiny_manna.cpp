#include "params.h"

#include <array>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstring>
#include <immintrin.h> // AVX2

#if N < 65536u
typedef uint16_t MannaSizeType;
#define AVX_WIDTH 16u
#else
typedef uint32_t MannaSizeType;
#define AVX_WIDTH 8u
#endif

typedef MannaSizeType Manna_Array[N];

// global state for xorshift rng
static uint32_t xrng_state = 0xdeadbeef;
static uint32_t xrng_bit_buffer = 0;
static uint8_t  xrng_bit_count = 0;

void xrng_init(uint32_t seed) {
    xrng_state = seed;
    xrng_bit_count = 0;
}

// Generate a 32-bit random number using xorshift32
uint32_t xrng_next() {
    uint32_t x = xrng_state;
    x ^= x << 13u;
    x ^= x >> 17u;
    x ^= x << 5u;
    xrng_state = x;
    return x;
}

// Generate a single random bit from a 32-bit buffered value
inline bool xrng_next_bit() {
    #if N >= 16384
    if (xrng_bit_count == 0) {
        xrng_bit_buffer = xrng_next();
        xrng_bit_count = 32;
    }
    bool bit = xrng_bit_buffer & 1u;
    xrng_bit_buffer >>= 1;
    --xrng_bit_count;
    return bit;
    #else
    return xrng_next() & 1u;
    #endif
}

// Vectorized version: generate 4 random 32-bit numbers at once
void xrng_next_vector(__m128i& out) {
    uint32_t r[4];
    for (int i = 0; i < 4; ++i) {
        r[i] = xrng_next();
    }
    out = _mm_set_epi32(r[3], r[2], r[1], r[0]);
}

// CONDICION INICIAL ---------------------------------------------------------------
static void inicializacion(Manna_Array &__restrict__ h)
{
    for (MannaSizeType i = 0u; i < N; ++i)
    {
        h[i] = static_cast<MannaSizeType>((i + 1) * DENSITY) - static_cast<MannaSizeType>(i * DENSITY);
    }
}

#ifdef DEBUG
static void progreso(const Manna_Array &__restrict__ h, std::ostream &__restrict__ output_file = std::cout)
{
    uint32_t granos = 0u;
    uint32_t granos_activos = 0u;

    for (MannaSizeType i = 0u; i < N; ++i)
    {
        output_file << h[i] << " ";
        granos += h[i];
        granos_activos += (h[i] > 1u);
    }
    output_file << granos << " - " << granos_activos << "\n";
}
#endif

// CONDICION INICIAL ---------------------------------------------------------------
static void desestabilizacion_inicial(Manna_Array &__restrict__ h)
{
    std::vector<MannaSizeType> index_a_incrementar;
    for (MannaSizeType i = 0u; i < N; ++i)
    {
        if (h[i] == 1)
        {
            h[i] = 0u;
            MannaSizeType j = (i + 2u * xrng_next_bit() - 1u) & (N - 1u);
            index_a_incrementar.push_back(j);
        }
    }
    for (MannaSizeType i = 0u; i < index_a_incrementar.size(); ++i)
    {
        h[index_a_incrementar[i]] += 1u;
    }
}

// DESCARGA DE ACTIVOS Y UPDATE --------------------------------------------------------
static MannaSizeType descargar(Manna_Array &__restrict__ h, Manna_Array &__restrict__ dh, uint32_t &__restrict__ granos_procesados)
{
    std::memset(dh, 0, sizeof(Manna_Array));

    for (MannaSizeType i = 0u; i < N; ++i)
    {
        if (h[i] > 1)
        {
            for (uint16_t j = 0u; j < h[i]; ++j)
            {
                MannaSizeType k = (i + 2u * xrng_next_bit() - 1u + N) & (N - 1u);
                ++dh[k];
            }
            granos_procesados += h[i];
            h[i] = 0u;
        }
    }

    MannaSizeType nroactivos = 0u;
    for (MannaSizeType i = 0u; i < N; ++i)
    {
        h[i] += dh[i];
        nroactivos += (h[i] > 1u);
    }

    return nroactivos;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    xrng_init(SEED);  // Initialize RNG
    Manna_Array h, dh;
    MannaSizeType activity;
    uint32_t t = 0u;
    uint32_t granos_procesados = 0u;

    inicializacion(h);
    #ifdef DEBUG
    std::ofstream output_file("sand.dat");
    progreso(h, output_file);
    #endif

    desestabilizacion_inicial(h);
    #ifdef DEBUG
    progreso(h, output_file);
    #endif

    do {
        activity = descargar(h, dh, granos_procesados);
        #ifdef DEBUG
        progreso(h, output_file);
        #endif
        ++t;
    } while (activity > 0u && t < NSTEPS);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Tiempo de procesamiento (s): " << static_cast<double>(duration.count()) / 1e6 << "\n";
    std::cout << "Granos procesados: " << granos_procesados << "\n";
    std::cout << "Granos/us: " << static_cast<double>(granos_procesados) / duration.count() << "\n";
    return 0;
}

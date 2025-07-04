#include "params.h"

#include <array>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include <cstring>

#if N < 65536
typedef uint16_t MannaSizeType;
#else
typedef uint32_t MannaSizeType;
#endif

static uint32_t rng_state = SEED;
typedef MannaSizeType Manna_Array[N];

#if N < 16384
static inline bool fast_rand()
{
    uint32_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rng_state = x;
    return x & 1;
}
#else
static uint32_t current_random_bits = 0u;
static uint32_t bits_used = 32u;
static inline bool fast_rand()
{
    if (bits_used >= 32)
    {
        uint32_t x = rng_state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        current_random_bits = rng_state = x;
        bits_used = 0;
    }

    bool bit = (current_random_bits >> bits_used) & 1;
    bits_used++;

    return bit;
}
#endif

// CONDICION INICIAL ---------------------------------------------------------------
static void inicializacion(Manna_Array &__restrict__ h)
{
    for (MannaSizeType i = 0u; i < N; ++i)
    {
        h[i] = (MannaSizeType)((i + 1) * DENSITY) - (MannaSizeType)(i * DENSITY);
    }
}

#ifdef DEBUG
static void progreso(const Manna_Array &h, std::ostream &output_file = std::cout)
{
    uint granos = 0;
    uint granos_activos = 0;

    for (int i = 0; i < N; ++i)
    {
        output_file << h[i] << " ";
        granos += h[i];
        granos_activos += (h[i] > 1);
    }
    output_file << "(" << granos << " " << granos_activos << " " << granos * 1.0 / N << " " << DENSITY << ")\n";
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
            h[i] = 0;
            int j = (i + 2 * fast_rand() - 1) & (N - 1);
            index_a_incrementar.push_back(j);
        }
    }
    for (MannaSizeType i = 0u; i < index_a_incrementar.size(); ++i)
    {
        h[index_a_incrementar[i]] += 1;
    }
}

// DESCARGA DE ACTIVOS Y UPDATE --------------------------------------------------------
static MannaSizeType descargar(Manna_Array &__restrict__ h, Manna_Array &__restrict__ dh, uint32_t &__restrict__ granos_procesados)
{
    memset(dh, 0, N * sizeof(MannaSizeType));

    for (MannaSizeType i = 0u; i < N; ++i)
    {
        if (h[i] > 1)
        {
            for (uint16_t j = 0u; j < h[i]; ++j)
            {
                MannaSizeType k = (i + 2 * fast_rand() - 1 + N) & (N - 1);
                ++dh[k];
            }
            granos_procesados += h[i];
            h[i] = 0u;
        }
    }

    MannaSizeType nroactivos = 0u;
    for (MannaSizeType i = 0; i < N; ++i)
    {
        h[i] += dh[i];
        nroactivos += (h[i] > 1);
    }

    return nroactivos;
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    Manna_Array h, dh;
    MannaSizeType activity;
    uint32_t t = 0u;
    uint32_t granos_procesados = 0u;

    inicializacion(h);
    desestabilizacion_inicial(h);

    do
    {
        activity = descargar(h, dh, granos_procesados);
        ++t;
    } while (activity > 0 && t < NSTEPS);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Tiempo de procesamiento (s): " << duration.count() / 1e6 << "\n";
    std::cout << "Granos procesados: " << granos_procesados << "\n";
    std::cout << "Granos/us: " << static_cast<double>(granos_procesados) / duration.count() << "\n";
    return 0;
}

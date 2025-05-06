#include "params.h"
#include "xorshift.hpp"

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <immintrin.h>

#if N < 65536u
typedef uint16_t MannaSizeType;
#else
typedef uint32_t MannaSizeType;
#endif

#define VEC_WIDTH 32u
typedef uint8_t MannaItemType;
typedef MannaItemType Manna_Array[N];


// AVX2 CONSTANTS ----------------------------------------------------------
// ones
static const __m256i one = _mm256_set1_epi8(1);
static const __m256i v8  = _mm256_set1_epi8(8);
static const __m256i all_ones = _mm256_set1_epi8((char)0xFF);
static const __m256i zeros = _mm256_setzero_si256();
// masks
static inline __m256i gt1(const __m256i &__restrict__ h)
{
    return _mm256_cmpgt_epi8(h, one);
}

template <typename T = uint8_t>
void print_vector(__m256i vec) {
    constexpr size_t num_elems = sizeof(vec) / sizeof(T);
    T* data = reinterpret_cast<T*>(&vec);
    for (size_t i = 0; i < num_elems; ++i) {
        std::cout << static_cast<int>(data[i]) << " ";
    }
    std::cout << std::endl;
}

// print bit per bit of 256i vector
void print_bits(__m256i vec) {
    constexpr size_t num_elems = sizeof(vec) * 8;
    uint8_t* data = reinterpret_cast<uint8_t*>(&vec);
    for (size_t i = 0; i < num_elems; ++i) {
        std::cout << ((data[i / 8] >> (7 - (i % 8))) & 1);
    }
    std::cout << std::endl;
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
    uint16_t max = 0u;

    for (MannaSizeType i = 0u; i < N; ++i)
    {
        output_file << static_cast<int>(h[i]) << " ";
        if (h[i] > max)
        {
            max = h[i];
        }
        granos += h[i];
        granos_activos += (h[i] > 1u);
    }
    output_file << " g=" << granos << " act=" << granos_activos << " max=" << static_cast<int>(max) << "\n";
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
            MannaSizeType j = (i + 2u * xrng_next_bit_idx(0) - 1u) & (N - 1u);
            index_a_incrementar.push_back(j);
        }
    }
    for (MannaSizeType i = 0u; i < index_a_incrementar.size(); ++i)
    {
        h[index_a_incrementar[i]] += 1u;
    }
}


static MannaSizeType descargar(Manna_Array &__restrict__ h,
                               Manna_Array &__restrict__ lh,
                               Manna_Array &__restrict__ rh)
{
    // 1) Inicializo acumuladores
    std::memset(lh, 0, sizeof(Manna_Array));
    std::memset(rh, 0, sizeof(Manna_Array));

    // 2) Lookup‑tables en 128‑bit (dos copias en cada YMM)
    static const __m128i mask_lo_lut = _mm_setr_epi8(
        // bits por elemento 0..15: (1<<i)-1 saturando a 0xFF si i>=8
        0x00,0x01,0x03,0x07,0x0F,0x1F,0x3F,0x7F,
        (char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF
    );
    static const __m128i mask_hi_lut = _mm_setr_epi8(
        // bits extra (h-8) 0..15: (1<<i)-1 saturado a 0xFF si i>=8
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x01,0x03,0x07,0x0F,0x1F,0x3F,0x7F,(char)0xFF
    );
    static const __m128i pop4_lut    = _mm_setr_epi8(
        // popcount de nibble 0..15
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4
    );

    // Broadcast a YMM
    __m256i m_lo_tbl  = _mm256_broadcastsi128_si256(mask_lo_lut);
    __m256i m_hi_tbl  = _mm256_broadcastsi128_si256(mask_hi_lut);
    __m256i pop4_tbl  = _mm256_broadcastsi128_si256(pop4_lut);

    // 2) Vectorizado de 32 elementos de 8 bits
    for (MannaSizeType i = 0; i < N; i += VEC_WIDTH) {
        // a) Carga alturas h[i..i+31]
        __m256i h_vec = _mm256_loadu_si256((__m256i*)&h[i]);

        // b) Mascara de sites activos (h>1)
        __m256i active = _mm256_cmpgt_epi8(h_vec, one);

        // c) RNG bajo y alto: 16 bits para cada site
        __m256i rnd_lo = xrng256_next();
        __m256i rnd_hi = xrng256_next();

        // d) indice hi = saturating_sub(h,8)
        __m256i hi_idx = _mm256_subs_epu8(h_vec, v8);

        // e) Lookup mascaras low/high bits
        __m256i mask_lo = _mm256_shuffle_epi8(m_lo_tbl, h_vec);
        __m256i mask_hi = _mm256_shuffle_epi8(m_hi_tbl, hi_idx);

        // f) mascara bits usados
        __m256i used_lo = _mm256_and_si256(rnd_lo, mask_lo);
        __m256i used_hi = _mm256_and_si256(rnd_hi, mask_hi);

        // g) Popcount nibble-wise para low y hi
        auto pop_byte = [&](const __m256i &x) {
            __m256i lo = _mm256_and_si256(x, _mm256_set1_epi8(0x0F));
            __m256i hi = _mm256_and_si256(_mm256_srli_epi16(x,4), _mm256_set1_epi8(0x0F));
            __m256i cnt_lo = _mm256_shuffle_epi8(pop4_tbl, lo);
            __m256i cnt_hi = _mm256_shuffle_epi8(pop4_tbl, hi);
            return _mm256_add_epi8(cnt_lo, cnt_hi);
        };
        __m256i pop_lo = pop_byte(used_lo);
        __m256i pop_hi = pop_byte(used_hi);

        // h) Conteo total y mascara activa
        __m256i left_cnt  = _mm256_add_epi8(pop_lo, pop_hi);
        left_cnt = _mm256_and_si256(left_cnt, active);

        // i) Calculo right = h - left
        __m256i raw_r = _mm256_sub_epi8(h_vec, left_cnt);
        __m256i right_cnt = _mm256_and_si256(raw_r, active);

        // j) Almaceno en buffers uint8_t
        _mm256_storeu_si256((__m256i*)&lh[i], left_cnt);
        _mm256_storeu_si256((__m256i*)&rh[i], right_cnt);
    }

    // 3) Condiciones de frontera periodicas
    {
        uint8_t last_r = rh[N-1];
        std::memmove(rh+1, rh, (N-1)*sizeof(MannaItemType)); rh[0]=last_r;
        uint8_t first_l= lh[0];
        std::memmove(lh, lh+1, (N-1)*sizeof(MannaItemType)); lh[N-1]=first_l;
    }

    // 4) Reconstruccion y cuenta actividad
    MannaSizeType activity=0;
    for (MannaSizeType i=0;i<N;++i) {
        h[i] = lh[i] + rh[i] + (h[i]==1u ? 1u:0u);
        activity += (h[i]>1u);
    }
    return activity;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    xrng256_init(SEED);
    alignas(32) Manna_Array h, lh, rh;
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
        activity = descargar(h, lh, rh);
        #ifdef DEBUG
        progreso(h, output_file);
        #endif
        ++t;
    } while (t < NSTEPS && activity > 0u);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Steps taken: " << t << "\n";
    std::cout << "Tiempo de procesamiento (s): " << static_cast<double>(duration.count()) / 1e6 << "\n";
    std::cout << "Granos procesados: " << granos_procesados << "\n";
    std::cout << "Granos/us: " << static_cast<double>(granos_procesados) / duration.count() << "\n";
    return 0;
}
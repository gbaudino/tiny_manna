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
#define VEC_ITERATIONS (N / VEC_WIDTH)
typedef uint8_t MannaItemType;
typedef MannaItemType Manna_Array[N];


// AVX2 CONSTANTS ----------------------------------------------------------
static const __m256i one = _mm256_set1_epi8(1);
static const __m256i v8  = _mm256_set1_epi8(8);
static const __m256i all_ones = _mm256_set1_epi8((char)0xFF);

template <typename T = uint8_t>
void print_vector(__m256i vec) {
    constexpr size_t num_elems = sizeof(vec) / sizeof(T);
    T* data = reinterpret_cast<T*>(&vec);
    for (size_t i = 0; i < num_elems; ++i) {
        std::cout << static_cast<int>(data[i]) << " ";
    }
    std::cout << std::endl;
}

// AVX2 LUTs (broadcasted to 256-bit)
alignas(32) static const __m256i lo_tbl256 = _mm256_broadcastsi128_si256(
    _mm_setr_epi8(
        0x00,0x01,0x03,0x07,0x0F,0x1F,0x3F,0x7F,
        (char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF,
        (char)0xFF,(char)0xFF,(char)0xFF,(char)0xFF
    )
);
alignas(32) static const __m256i hi_tbl256 = _mm256_broadcastsi128_si256(
    _mm_setr_epi8(
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x01,0x03,0x07,0x0F,0x1F,0x3F,0x7F,(char)0xFF
    )
);
alignas(32) static const __m256i pop4_tbl256 = _mm256_broadcastsi128_si256(
    _mm_setr_epi8( 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 )
);

// Inline popcount for each byte in 256-bit vector
static inline __m256i popcount_byte_avx2(__m256i x) {
    __m256i lo = _mm256_and_si256(x, _mm256_set1_epi8(0x0F));
    __m256i hi = _mm256_and_si256(_mm256_srli_epi16(x, 4), _mm256_set1_epi8(0x0F));
    __m256i cnt_lo = _mm256_shuffle_epi8(pop4_tbl256, lo);
    __m256i cnt_hi = _mm256_shuffle_epi8(pop4_tbl256, hi);
    return _mm256_add_epi8(cnt_lo, cnt_hi);
}

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
                               Manna_Array &__restrict__ rh,
                               uint32_t &__restrict__ processed)
{
    // 1) Inicializo acumuladores
    std::memset(lh, 0, sizeof(Manna_Array));
    std::memset(rh, 0, sizeof(Manna_Array));

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
        __m256i mask_lo = _mm256_shuffle_epi8(lo_tbl256, h_vec);
        __m256i mask_hi = _mm256_shuffle_epi8(hi_tbl256, hi_idx);

        // f) mascara bits usados
        __m256i used_lo = _mm256_and_si256(rnd_lo, mask_lo);
        __m256i used_hi = _mm256_and_si256(rnd_hi, mask_hi);

        // g) Popcount nibble-wise para low y hi
        __m256i pop_lo = popcount_byte_avx2(used_lo);
        __m256i pop_hi = popcount_byte_avx2(used_hi);

        // h) Conteo total y mascara activa
        __m256i left_cnt  = _mm256_add_epi8(pop_lo, pop_hi);
        left_cnt = _mm256_and_si256(left_cnt, active);

        // i) Calculo right = h - left
        __m256i raw_r = _mm256_sub_epi8(h_vec, left_cnt);
        __m256i right_cnt = _mm256_and_si256(raw_r, active);

        // j) Almaceno en buffers uint8_t
        _mm256_storeu_si256((__m256i*)&lh[i], left_cnt);
        _mm256_storeu_si256((__m256i*)&rh[i], right_cnt);

        // k) Acumulo granos procesados
        __m256i pop_vec = _mm256_add_epi8(left_cnt, right_cnt);
        __m256i sums = _mm256_sad_epu8(pop_vec, _mm256_setzero_si256());
        uint64_t chunk = _mm256_extract_epi64(sums, 0)
                       + _mm256_extract_epi64(sums, 1)
                       + _mm256_extract_epi64(sums, 2)
                       + _mm256_extract_epi64(sums, 3);
        processed += chunk;
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
    uint32_t processed = 0u;

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
        activity = descargar(h, lh, rh, processed);
        #ifdef DEBUG
        progreso(h, output_file);
        #endif
        ++t;
    } while (t < NSTEPS && activity > 0u);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Steps taken: " << t << "\n";
    std::cout << "Tiempo de procesamiento (s): " << static_cast<double>(duration.count()) / 1e6 << "\n";
    std::cout << "Granos procesados: " << processed << "\n";
    std::cout << "Granos/us: " << static_cast<double>(processed) / duration.count() << "\n";
    return 0;
}
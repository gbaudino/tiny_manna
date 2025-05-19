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
#include <omp.h>

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

// CONDICION INICIAL ---------------------------------------------------------------
static void inicializacion(Manna_Array &__restrict__ h) {
    #pragma omp parallel for schedule(static)
    for (MannaSizeType i = 0u; i < N; ++i) {
        h[i] = static_cast<MannaSizeType>((i + 1) * DENSITY) - static_cast<MannaSizeType>(i * DENSITY);
    }
}

#ifdef DEBUG
static void progreso(const Manna_Array &__restrict__ h, std::ostream &__restrict__ output_file = std::cout) {
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

static void desestabilizacion_inicial(Manna_Array &__restrict__ h) {
    std::vector<MannaSizeType> index_a_incrementar;
    for (MannaSizeType i = 0u; i < N; ++i) {
        if (h[i] == 1) {
            h[i] = 0u;
            MannaSizeType j = (i + 2u * xrng_next_bit() - 1u) & (N - 1u);
            index_a_incrementar.push_back(j);
        }
    }
    for (MannaSizeType i = 0u; i < index_a_incrementar.size(); ++i) {
        h[index_a_incrementar[i]] += 1u;
    }
}


// Descarga paralelizada
// static bool descargar(Manna_Array &h,
//                       Manna_Array &lh,
//                       Manna_Array &rh,
//                       uint32_t  &processed)
// {
//     // buffers y contadores
//     std::memset(lh,0,sizeof(h));
//     std::memset(rh,0,sizeof(h));

//     // flag de actividad global y suma de granos
//     bool any_active = false;

//     #pragma omp parallel
//     {
//         // acumulador de granos por hilo
//         __m256i v_acc = _mm256_setzero_si256();
//         bool local_active = false;

//         // cada hilo procesa un tramo
//         #pragma omp for schedule(static) nowait
//         for (MannaSizeType i = 0; i < N; i += VEC_WIDTH) {
//             // 1) carga
//             __m256i h_vec = _mm256_loadu_si256((__m256i*)(h + i));
//             __m256i active = _mm256_cmpgt_epi8(h_vec, one);
//             // RNG
//             __m256i rnd_lo = xrng256_next();
//             __m256i rnd_hi = xrng256_next();
//             // vectores
//             __m256i hi_idx  = _mm256_subs_epu8(h_vec, v8);
//             __m256i mask_lo = _mm256_shuffle_epi8(lo_tbl256, h_vec);
//             __m256i mask_hi = _mm256_shuffle_epi8(hi_tbl256, hi_idx);
//             __m256i used_lo = _mm256_and_si256(rnd_lo, mask_lo);
//             __m256i used_hi = _mm256_and_si256(rnd_hi, mask_hi);
//             __m256i pop_lo  = popcount_byte_avx2(used_lo);
//             __m256i pop_hi  = popcount_byte_avx2(used_hi);
//             __m256i left_cnt  = _mm256_and_si256(_mm256_add_epi8(pop_lo,pop_hi), active);
//             __m256i right_cnt = _mm256_and_si256(_mm256_sub_epi8(h_vec,left_cnt), active);
//             // guarda
//             _mm256_storeu_si256((__m256i*)(lh + i), left_cnt);
//             _mm256_storeu_si256((__m256i*)(rh + i), right_cnt);
//             // acumula
//             __m256i tot = _mm256_add_epi8(left_cnt, right_cnt);
//             __m256i sums = _mm256_sad_epu8(tot, _mm256_setzero_si256());
//             v_acc = _mm256_add_epi64(v_acc, sums);
//             // chequear actividad
//             if (!_mm256_testz_si256(active, active)) local_active = true;
//         }

//         // reducc de granos procesados
//         uint64_t buf[4];
//         _mm256_storeu_si256((__m256i*)buf, v_acc);
//         uint32_t part = uint32_t(buf[0]+buf[1]+buf[2]+buf[3]);
//         #pragma omp atomic
//         processed += part;

//         // actualizar flag global
//         #pragma omp critical
//         any_active |= local_active;

//         // barrera para asegurar que todos terminaron
//         #pragma omp barrier

//         // sección única para fronteras periódicas
//         #pragma omp single
//         {
//             uint8_t last_r = rh[N-1];
//             std::memmove(rh+1, rh, (N-1)*sizeof(MannaItemType)); rh[0] = last_r;
//             uint8_t first_l = lh[0];
//             std::memmove(lh, lh+1, (N-1)*sizeof(MannaItemType)); lh[N-1] = first_l;
//         }

//         // reconstrucción paralela
//         #pragma omp for schedule(static)
//         for (MannaSizeType i = 0; i < N; ++i) {
//             h[i] = lh[i] + rh[i] + (h[i]==1u ? 1u : 0u);
//         }
//     }

//     return any_active;
// }

static bool descargar(Manna_Array &h,
                      Manna_Array &lh,
                      Manna_Array &rh,
                      uint32_t  &processed)
{
    // reset buffers
    std::memset(lh,0,sizeof(h));
    std::memset(rh,0,sizeof(h));

    bool any_active = false;
    uint32_t local_processed = 0u;

    // 1) vector loop paralelo con reducciones
    #pragma omp parallel for reduction(||:any_active) reduction(+:local_processed) schedule(static)
    for (MannaSizeType i = 0; i < N; i += VEC_WIDTH) {
        __m256i h_vec    = _mm256_loadu_si256((__m256i*)(h + i));
        __m256i active   = _mm256_cmpgt_epi8(h_vec, one);
        // RNG local por hilo
        __m256i rnd_lo   = xrng256_next();
        __m256i rnd_hi   = xrng256_next();
        __m256i hi_idx   = _mm256_subs_epu8(h_vec, v8);
        __m256i mask_lo  = _mm256_shuffle_epi8(lo_tbl256, h_vec);
        __m256i mask_hi  = _mm256_shuffle_epi8(hi_tbl256, hi_idx);
        __m256i used_lo  = _mm256_and_si256(rnd_lo, mask_lo);
        __m256i used_hi  = _mm256_and_si256(rnd_hi, mask_hi);
        __m256i pop_lo   = popcount_byte_avx2(used_lo);
        __m256i pop_hi   = popcount_byte_avx2(used_hi);
        __m256i left_cnt  = _mm256_and_si256(_mm256_add_epi8(pop_lo,pop_hi), active);
        __m256i right_cnt = _mm256_and_si256(_mm256_sub_epi8(h_vec,left_cnt), active);
        _mm256_storeu_si256((__m256i*)(lh + i), left_cnt);
        _mm256_storeu_si256((__m256i*)(rh + i), right_cnt);
        // actualizar reducciones
        uint64_t buf[4];
        __m256i tot = _mm256_add_epi8(left_cnt, right_cnt);
        __m256i sums = _mm256_sad_epu8(tot, _mm256_setzero_si256());
        _mm256_storeu_si256((__m256i*)buf, sums);
        local_processed += uint32_t(buf[0]+buf[1]+buf[2]+buf[3]);
        any_active |= !_mm256_testz_si256(active, active);
    }
    // sumar al contador global
    processed += local_processed;

    // 2) frontera periódica (único hilo)
    #pragma omp single
    {
        uint8_t last_r = rh[N-1];
        std::memmove(rh+1, rh, (N-1)*sizeof(MannaItemType)); rh[0] = last_r;
        uint8_t first_l = lh[0];
        std::memmove(lh, lh+1, (N-1)*sizeof(MannaItemType)); lh[N-1] = first_l;
    }

    // 3) reconstrucción paralela
    for (MannaSizeType i = 0; i < N; ++i) {
        h[i] = lh[i] + rh[i] + (h[i]==1u ? 1u:0u);
    }

    return any_active;
}

int main() {
    omp_set_num_threads(THREADS);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        xrng256_thread_init(SEED, tid);
    }

    auto start = std::chrono::high_resolution_clock::now();
    
    alignas(32) Manna_Array h, lh, rh;
    bool activity;
    uint32_t t = 0u;
    uint32_t processed_global = 0u;

    inicializacion(h);
    #ifdef DEBUG
    std::ofstream output_file("sand.dat");
    progreso(h, output_file);
    #endif

    desestabilizacion_inicial(h);
    #ifdef DEBUG
    progreso(h, output_file);
    #endif

    bool any_active_global = false;

    #pragma omp parallel shared(h, processed_global, any_active_global, t)
    {
        int tid = omp_get_thread_num();
        size_t chunk = N / THREADS;
        size_t start_idx = tid * chunk;
        size_t end_idx = start_idx + chunk;

        // buffers locales alineados
        std::vector<MannaItemType> lh_loc(chunk), rh_loc(chunk);

        bool any_active_local;
        uint32_t processed_local;

        do {
            // reset locales y globales (por hilo)
            any_active_local = false;
            processed_local = 0;

            // 1) Loop vectorizado sobre el chunk
            for (size_t off = 0; off < chunk; off += VEC_WIDTH) {
                size_t i = start_idx + off;
                __m256i h_vec  = _mm256_loadu_si256((__m256i*)(h + i));
                __m256i active = _mm256_cmpgt_epi8(h_vec, one);

                __m256i rnd_lo = xrng256_next();
                __m256i rnd_hi = xrng256_next();
                __m256i hi_idx = _mm256_subs_epu8(h_vec, v8);
                __m256i mask_lo = _mm256_shuffle_epi8(lo_tbl256, h_vec);
                __m256i mask_hi = _mm256_shuffle_epi8(hi_tbl256, hi_idx);
                __m256i used_lo = _mm256_and_si256(rnd_lo, mask_lo);
                __m256i used_hi = _mm256_and_si256(rnd_hi, mask_hi);
                __m256i pop_lo = popcount_byte_avx2(used_lo);
                __m256i pop_hi = popcount_byte_avx2(used_hi);
                __m256i left_cnt  = _mm256_and_si256(_mm256_add_epi8(pop_lo, pop_hi), active);
                __m256i right_cnt = _mm256_and_si256(_mm256_sub_epi8(h_vec, left_cnt), active);

                _mm256_storeu_si256((__m256i*)(lh_loc.data() + off), left_cnt);
                _mm256_storeu_si256((__m256i*)(rh_loc.data() + off), right_cnt);

                // conteo
                uint64_t buf[4];
                __m256i tot  = _mm256_add_epi8(left_cnt, right_cnt);
                __m256i sums = _mm256_sad_epu8(tot, _mm256_setzero_si256());
                _mm256_storeu_si256((__m256i*)buf, sums);
                processed_local += uint32_t(buf[0] + buf[1] + buf[2] + buf[3]);
                any_active_local |= !_mm256_testz_si256(active, active);
            }

            // 2) Periodicidad dentro del chunk (memmove)
            uint8_t last_r = rh_loc[chunk - 1];
            std::memmove(rh_loc.data() + 1, rh_loc.data(), (chunk - 1) * sizeof(MannaItemType));
            rh_loc[0] = last_r;
            uint8_t first_l = lh_loc[0];
            std::memmove(lh_loc.data(), lh_loc.data() + 1, (chunk - 1) * sizeof(MannaItemType));
            lh_loc[chunk - 1] = first_l;

            // 3) Reconstrucción global por bloques
            for (size_t off = 0; off < chunk; ++off) {
                size_t i = start_idx + off;
                h[i] = lh_loc[off] + rh_loc[off] + (h[i] == 1u ? 1u : 0u);
            }

            // 4) Reducciones atómicas/implicitas
            #pragma omp critical
            {
                processed_global += processed_local;
                any_active_global |= any_active_local;
            }

            

            #pragma omp single
            {
                #ifdef DEBUG
                progreso(h, output_file);
                #endif
                ++t;
            }


        } while (any_active_global && t < NSTEPS);
    }



    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Pasos: " << t << "\n";
    std::cout << "Tiempo de procesamiento (s): " << static_cast<double>(duration.count()) / 1e6 << "\n";
    std::cout << "Granos procesados: " << processed_global << "\n";
    std::cout << "Granos/us: " << static_cast<double>(processed_global) / duration.count() << "\n";
    return 0;
}
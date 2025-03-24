#include "params.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>

typedef std::array<int, N> Manna_Array;


// CONDICION INICIAL ---------------------------------------------------------------
static void inicializacion(Manna_Array& h)
{
    for (int i = 0; i < N; ++i) {
        h[i] = static_cast<int>((i + 1) * DENSITY) - static_cast<int>(i * DENSITY);
    }
}


#ifdef DEBUG
static void progreso(const Manna_Array& h, std::ostream& output_file = std::cout)
{
    int granos = 0;
    int granos_activos = 0;

    for (int i = 0; i < N; ++i) {
        output_file << h[i] << " ";
        granos += h[i];
        granos_activos += (h[i] > 1);
    }
    output_file << "(" <<  granos << " " << granos_activos << " " << granos * 1.0 / N << " " << DENSITY << ")\n";
}
#endif


// CONDICION INICIAL ---------------------------------------------------------------
static void desestabilizacion_inicial(Manna_Array& h)
{
    std::vector<int> index_a_incrementar;
    for (int i = 0; i < N; ++i) {
        if (h[i] == 1) {
            h[i] = 0;
            int j = i + 2 * (rand() % 2) - 1; // izquierda o derecha

            if (j == N) {
                j = 0;
            } else if (j == -1) {
                j = N - 1;
            }

            index_a_incrementar.push_back(j);
        }
    }
    for (unsigned int i = 0; i < index_a_incrementar.size(); ++i) {
        h[index_a_incrementar[i]] += 1;
    }
}


// DESCARGA DE ACTIVOS Y UPDATE --------------------------------------------------------
static unsigned int descargar(Manna_Array& h, Manna_Array& dh)
{
    dh.fill(0);

    for (int i = 0; i < N; ++i) {
        if (h[i] > 1) {
            for (int j = 0; j < h[i]; ++j) {
                int k = (i + 2 * (rand() % 2) - 1 + N) % N;
                ++dh[k];
            }
            h[i] = 0;
        }
    }

    unsigned int nroactivos = 0;
    for (int i = 0; i < N; ++i) {
        h[i] += dh[i];
        nroactivos += (h[i] > 1);
    }

    return nroactivos;
}

int main()
{
    srand(SEED);
    Manna_Array h, dh;
    int activity;
    int t = 0;
    unsigned long long granos_procesados = 0;
    
    inicializacion(h);
    #ifdef DEBUG
    std::ofstream output_file("sand.dat");
    progreso(h, output_file);
    #endif

    auto start = std::chrono::high_resolution_clock::now();

    desestabilizacion_inicial(h);
    #ifdef DEBUG
    progreso(h, output_file);
    #endif

    do {
        activity = descargar(h, dh);
        #ifdef DEBUG
        progreso(h, output_file);
        #endif
        granos_procesados += std::accumulate(h.begin(), h.end(), 0ULL);
        ++t;
    } while (activity > 0 && t < NSTEPS);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double granos_us = static_cast<double>(granos_procesados) / elapsed_us.count();

    #ifdef DEBUG
    auto elapsed_s = std::chrono::duration<double>(end - start);
    std::cout << "Tiempo total: " << elapsed_s.count() << " segundos\n";
    std::cout << "Granos: " << granos_procesados/1E6 << "M\n";
    std::cout << (t == NSTEPS ? "Se terminó el tiempo" : "La actividad decayó a cero") << "\n";
    #endif

    std::cout << "Granos/µs: " << granos_us << "\n";
    return 0;
}

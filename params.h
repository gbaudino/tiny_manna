#pragma once

#include <ctime>

#ifndef N
constexpr size_t N = 4096;
#endif

#ifndef DENSITY
#define DENSITY 0.8924
#endif

#ifndef NSTEPS
#define NSTEPS 1000000u
#endif

#ifndef SEED
#define SEED 1 // 1 if you want fixed seed
#endif

#ifndef THREADS
#define THREADS 1
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 256
#endif

// #ifndef DEBUG
// #define DEBUG 1
// #endif

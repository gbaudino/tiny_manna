#pragma once

#include <ctime>

// number of sites
#ifndef N
#define N 16384
#endif

// number of sites
#ifndef DENSITY
#define DENSITY 0.8924
#endif

// number of temporal steps
#ifndef NSTEPS
#define NSTEPS 1000000
#endif

#ifndef SEED
#define SEED 1 //time(nullptr) // 1 if you want fixed seed
#endif

// #ifndef DEBUG
// #define DEBUG 1
// #endif

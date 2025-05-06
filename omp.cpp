// openmp_introduction.cpp
// Introduction to using OpenMP in C++
// Compile with: g++ -std=c++17 -fopenmp openmp_introduction.cpp -o openmp_example
// Run with: OMP_NUM_THREADS=<num_threads> ./openmp_example

#include <iostream>
#include <omp.h>
#include <vector>

// Utility: print header with thread info
void print_header(const std::string &title) {
    #pragma omp single
    {
        std::cout << "===== " << title << " =====" << std::endl;
    }
}

int main() {
    // 1. Basic Parallel Region
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        #pragma omp critical
        {
            std::cout << "Hello from thread " << tid << " of " << nthreads << std::endl;
        }
    }
    std::cout << std::endl;

    // 2. Parallel For with default scheduling
    const int N = 10;
    std::vector<int> v(N);

    print_header("Parallel For");
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        v[i] = i * i;
        #pragma omp critical
        std::cout << "Thread " << omp_get_thread_num()
                  << " computed v[" << i << "] = " << v[i] << std::endl;
    }
    std::cout << std::endl;

    // 3. Reduction Example (sum of v)
    print_header("Reduction");
    long sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; ++i) {
        sum += v[i];
    }
    std::cout << "Sum of squares = " << sum << std::endl << std::endl;

    // 4. Scheduling clauses
    print_header("Scheduling Clauses");
    #pragma omp parallel for schedule(static, 2)
    for (int i = 0; i < N; ++i) {
        #pragma omp critical
        std::cout << "static schedule: thread " << omp_get_thread_num() << " -> i=" << i << std::endl;
    }
    std::cout << std::endl;
    #pragma omp parallel for schedule(dynamic, 3)
    for (int i = 0; i < N; ++i) {
        #pragma omp critical
        std::cout << "dynamic schedule: thread " << omp_get_thread_num() << " -> i=" << i << std::endl;
    }
    std::cout << std::endl;

    // 5. Sections Example
    print_header("Sections");
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            std::cout << "Section 1 executed by thread " << omp_get_thread_num() << std::endl;
        }
        #pragma omp section
        {
            std::cout << "Section 2 executed by thread " << omp_get_thread_num() << std::endl;
        }
        #pragma omp section
        {
            std::cout << "Section 3 executed by thread " << omp_get_thread_num() << std::endl;
        }
    }
    std::cout << std::endl;

    // 6. Single and Master
    print_header("Single vs Master");
    #pragma omp parallel
    {
        #pragma omp single
        {
            std::cout << "Only one thread (not necessarily master) executes this: "
                      << omp_get_thread_num() << std::endl;
        }

        #pragma omp master
        {
            std::cout << "Only master thread (tid=0) executes this." << std::endl;
        }
    }
    std::cout << std::endl;

    // 7. Atomic Update
    print_header("Atomic");
    long counter = 0;
    #pragma omp parallel for
    for (int i = 0; i < 10000; ++i) {
        #pragma omp atomic
        ++counter;
    }
    std::cout << "Counter = " << counter << std::endl << std::endl;

    // 8. Barrier and ordered
    print_header("Barrier and Ordered");
    #pragma omp parallel for ordered
    for (int i = 0; i < 8; ++i) {
        #pragma omp ordered
        {
            std::cout << "Ordered output: " << i << " by thread "
                      << omp_get_thread_num() << std::endl;
        }
    }
    std::cout << std::endl;

    // 9. Tasks (OpenMP 3.0+)
    print_header("Tasks");
    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 1; i <= 5; ++i) {
                #pragma omp task firstprivate(i)
                {
                    std::cout << "Task " << i << " executed by thread "
                              << omp_get_thread_num() << std::endl;
                }
            }
        }
    }
    std::cout << std::endl;

    // 10. Threadprivate data
    static int thread_data = 0;
    #pragma omp threadprivate(thread_data)
    #pragma omp parallel
    {
        thread_data = omp_get_thread_num();
        #pragma omp critical
        std::cout << "Thread-private data for thread "
                  << omp_get_thread_num() << " = " << thread_data << std::endl;
    }

    return 0;
}

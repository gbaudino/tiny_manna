#!/bin/bash

# Check if all required arguments are provided
if [ $# -ne 4 ]; then
    echo "Usage: $0 <N> <RUNS> <CXX> <FLAGS>"
    echo "Example: $0 16384 2 g++ \"-O3 -march=native\""
    exit 1
fi

# Parse input arguments
N=$1
RUNS=$2
CXX=$3
FLAGS=$4

# Verify compiler exists
if ! command -v $CXX &> /dev/null; then
    echo "Compiler $CXX not found!"
    exit 1
fi

# Create results directory
RESULTS_DIR="benchmarks"
mkdir -p $RESULTS_DIR

# Modify params.h with the new N
sed -i "s/^#define N .*/#define N $N/" params.h

# Prepare results file
RESULTS_FILE="$RESULTS_DIR/${CXX}_${N}.csv"
echo "N,Compiler,Flags,Best_Grains/µs,Total_Time(s),Memory_KB" > $RESULTS_FILE

# Clean and compile
make clean > /dev/null 2>&1
OPTFLAGS="$FLAGS" CXX=$CXX make > /dev/null 2>&1

if [ ! -f tiny_manna ]; then
    echo "Compilation error with $CXX and flags $FLAGS"
    exit 1
fi

# Benchmark execution
echo "Running benchmark with:"
echo " - N: $N"
echo " - Runs: $RUNS"
echo " - Compiler: $CXX"
echo " - Flags: $FLAGS"

BEST_IPS=0
BEST_TIME=99999
BEST_MEM=999999999

for ((j=1; j<=RUNS; j++)); do
    echo -n "Run $j/$RUNS... "

    OUTPUT=$(perf stat -e instructions /usr/bin/time -v ./tiny_manna 2>&1)

    EXEC_TIME=$(echo "$OUTPUT" | grep "seconds time elapsed" | awk '{print $1}' | sed 's/,/./')
    INSTRUCTIONS=$(echo "$OUTPUT" | grep "instructions" | awk '{print $1}' | tr -d '.')
    MEM_USED=$(echo "$OUTPUT" | grep "Maximum resident set size" | awk '{print $6}')
    M_IPS=$(echo "scale=5; $INSTRUCTIONS / $EXEC_TIME / 1000000" | bc -l)
    
    echo "Instructions: $INSTRUCTIONS | IPS: $M_IPS M | Time: $EXEC_TIME s | Memory: ${MEM_USED} KB"
    
    if (( $(echo "$M_IPS > $BEST_IPS" | bc -l) )); then
        BEST_IPS=$M_IPS
    fi
    
    if (( $(echo "$EXEC_TIME < $BEST_TIME" | bc -l) )); then
        BEST_TIME=$EXEC_TIME
    fi
    
    if (( MEM_USED < BEST_MEM )); then
        BEST_MEM=$MEM_USED
    fi
done

# Write results to CSV
echo "$N,$CXX,\"$FLAGS\",$BEST_IPS,$BEST_TIME,$BEST_MEM" >> $RESULTS_FILE

# Display final results
echo -e "\nBest results:"
echo " - IPS: $BEST_IPS M"
echo " - Total Time (s): $BEST_TIME"
echo " - Memory (KB): $BEST_MEM"
echo -e "\nResults saved in $RESULTS_FILE"

# Clean up
make clean > /dev/null 2>&1
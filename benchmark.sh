#!/bin/bash
# Directorio para resultados
RESULTS_DIR="benchmarks"
mkdir -p $RESULTS_DIR

# Parámetros de entrada
DEVICE=${1:-"default"}        # Identificador de dispositivo (primer parámetro o "default" por defecto)
OPTIMIZATION=${2:-"std"}      # Tipo de optimización (segundo parámetro o "std" por defecto)

# Archivo de resultados
RESULTS_FILE="$RESULTS_DIR/${DEVICE}_${OPTIMIZATION}.csv"

# Definir array de valores N a probar y sus respectivas cantidades de ejecuciones
# Formato: N:RUNS (Ejemplo: 16384:5 significa probar N=16384 con 5 ejecuciones)
N_RUNS_ARRAY=(
  "512:8"
  "1024:4"
  "2048:4"
  "4096:4"
  "8192:2"
  "16384:1"
  #"32768:1"
)

# Lista de compiladores a probar
COMPILERS=(
  "g++"
  "clang++"
  #"icpx"
  #"nvc++"
)

# Flags de compilación a probar
FLAGS_LIST=(
  ""
  "-O2"
  "-O3 -march=native"
 #"-Ofast"
  "-O3 -funroll-loops"
)

# Encabezado del archivo CSV
echo "N,Compiler,Flags,Grains/µs,Total_Time(s),Memory_KB" > $RESULTS_FILE

for N_RUNS in "${N_RUNS_ARRAY[@]}"; do
  N=$(echo $N_RUNS | cut -d: -f1)
  RUNS=$(echo $N_RUNS | cut -d: -f2)
  
  echo "====================================="
  echo "Probando con N = $N ($RUNS ejecuciones)"
  echo "====================================="
  
  sed -i "s/^#define N .*/#define N $N/" params.h
  
  for CXX in "${COMPILERS[@]}"; do
    if ! command -v $CXX &> /dev/null; then
      echo "Compilador $CXX no encontrado, omitiendo..."
      continue
    fi
    
    for FLAGS in "${FLAGS_LIST[@]}"; do
      echo "Probando con $CXX y flags: $FLAGS"
      make clean > /dev/null 2>&1
      OPTFLAGS="$FLAGS" CXX=$CXX make > /dev/null 2>&1
      
      if [ ! -f tiny_manna ]; then
        echo "Error compilando con $CXX y flags $FLAGS, omitiendo..."
        continue
      fi
      
      echo "Ejecutando $RUNS veces..."
      BEST_GRAINS_PER_US=0
      BEST_TIME=99999
      BEST_MEM=999999999
      
      for ((j=1; j<=$RUNS; j++)); do
        echo -n "Ejecución $j/$RUNS... "
        
        OUTPUT=$(perf stat -e instructions /usr/bin/time -v ./tiny_manna 2>&1)
        
        EXEC_TIME=$(echo "$OUTPUT" | grep "seconds time elapsed" | awk '{print $1}' | sed 's/,/./')
        MEM_USED=$(echo "$OUTPUT" | grep "Maximum resident set size" | awk '{print $6}')
        GRAINS=$(echo "$OUTPUT" | grep "Granos activos:" | awk '{print $3}')
        GRAINS_PER_US=$(echo "scale=5; $GRAINS / ($EXEC_TIME * 1000000)" | bc -l)
        
        echo " Grains: $GRAINS | Grains/µs: $GRAINS_PER_US | Time: $EXEC_TIME s | Memory: ${MEM_USED} KB"
        
        if (( $(echo "$GRAINS_PER_US > $BEST_GRAINS_PER_US" | bc -l) )); then
          BEST_GRAINS_PER_US=$GRAINS_PER_US
        fi
        
        if (( $(echo "$EXEC_TIME < $BEST_TIME" | bc -l) )); then
          BEST_TIME=$EXEC_TIME
        fi
        
        if (( MEM_USED < BEST_MEM )); then
          BEST_MEM=$MEM_USED
        fi
      done
      
      echo "$N,$CXX,\"$FLAGS\",$BEST_GRAINS_PER_US,$BEST_TIME,$BEST_MEM" >> $RESULTS_FILE
      
      echo "Mejor resultado para N=$N, $CXX con flags $FLAGS:"
      echo " - Grains/µs: $BEST_GRAINS_PER_US"
      echo " - Tiempo total (s): $BEST_TIME"
      echo " - Memoria (KB): $BEST_MEM"
      echo "-----------------------------------------"
    done
  done
done

make clean > /dev/null 2>&1

echo -e "\nResultados guardados en $RESULTS_FILE"
cat $RESULTS_FILE

echo -e "\nTodos los resultados guardados en $RESULTS_FILE"

#!/bin/bash
# Directorio para resultados
RESULTS_DIR="benchmarks"
mkdir -p $RESULTS_DIR

# Parámetros de entrada
DEVICE=${1:-"default"}        # Identificador de dispositivo
OPTIMIZATION=${2:-"std"}      # Tipo de optimización

# Archivo de resultados
RESULTS_FILE="$RESULTS_DIR/${DEVICE}_${OPTIMIZATION}.csv"

# Definir array de valores N a probar y sus respectivas cantidades de ejecuciones
N_RUNS_ARRAY=(
  "4096:4"
  "8192:4"
  "16384:2"
  "32768:2"
  "65536:1"
  "131072:1"
  "262144:1"
)

# Lista de compiladores y flags combinados
COMPILER_FLAGS_LIST=(
  "g++: -O3 -march=native -funroll-loops"
  "clang++: -O2"
  "clang++: -O3 -march=native -funroll-loops -ffast-math"
  "icpx: -xHost -O3 -mavx2"
)

# Si el archivo no existe, crear el encabezado
if [ ! -f "$RESULTS_FILE" ]; then
  echo "N,Compiler,Flags,Grains/µs,Total_Time(s),Memory_KB" > "$RESULTS_FILE"
fi

for N_RUNS in "${N_RUNS_ARRAY[@]}"; do
  N=$(echo $N_RUNS | cut -d: -f1)
  RUNS=$(echo $N_RUNS | cut -d: -f2)

  echo "====================================="
  echo "Probando con N = $N ($RUNS ejecuciones)"
  echo "====================================="
  
  sed -i "s/^#define N .*/#define N $N/" params.h

  for COMPILER_FLAGS in "${COMPILER_FLAGS_LIST[@]}"; do
    CXX=$(echo $COMPILER_FLAGS | cut -d: -f1)
    FLAGS=$(echo $COMPILER_FLAGS | cut -d: -f2-)

    if ! command -v $CXX &> /dev/null; then
      echo "Compilador $CXX no encontrado, omitiendo..."
      continue
    fi

    # Verificar si la combinación ya existe en el archivo CSV
    if grep -q "^$N,$CXX,\"$FLAGS\"," "$RESULTS_FILE"; then
      echo "Resultados para N=$N, $CXX con flags \"$FLAGS\" ya existen. Omitiendo..."
      continue
    fi

    BEST_GRAINS_PER_US=0
    BEST_TIME=99999
    BEST_MEM=999999999

    for ((j=1; j<=$RUNS; j++)); do
      echo "Compilando con $CXX y flags: $FLAGS (Ejecución $j/$RUNS)"
      make clean > /dev/null 2>&1
      OPTFLAGS="$FLAGS" CXX=$CXX make > /dev/null 2>&1

      if [ ! -f tiny_manna ]; then
        echo "Error compilando con $CXX y flags $FLAGS, omitiendo..."
        continue
      fi

      echo -n "Ejecución $j/$RUNS... "

      OUTPUT=$(./tiny_manna 2>&1)

      EXEC_TIME=$(echo "$OUTPUT" | grep "Tiempo de procesamiento" | awk '{print $5}')
      GRAINS_PER_US=$(echo "$OUTPUT" | grep "Granos/us" | awk '{print $2}')
      MEM_USED=$(ps -o rss= -p $$ | tr -d ' ')

      echo " Grains/µs: $GRAINS_PER_US | Time: $EXEC_TIME s | Memory: ${MEM_USED} KB"

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

    # Guardar solo si no existe en el archivo
    echo "$N,$CXX,\"$FLAGS\",$BEST_GRAINS_PER_US,$BEST_TIME,$BEST_MEM" >> "$RESULTS_FILE"

    echo "Mejor resultado para N=$N, $CXX con flags $FLAGS:"
    echo " - Grains/µs: $BEST_GRAINS_PER_US"
    echo " - Tiempo total (s): $BEST_TIME"
    echo " - Memoria (KB): $BEST_MEM"
    echo "-----------------------------------------"
  done
done

make clean > /dev/null 2>&1

echo -e "\nResultados guardados en $RESULTS_FILE"
cat "$RESULTS_FILE"

echo -e "\nTodos los resultados guardados en $RESULTS_FILE"

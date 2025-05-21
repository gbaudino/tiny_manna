#!/bin/bash

# Directorio para resultados
RESULTS_DIR="benchmarks"
mkdir -p "$RESULTS_DIR"

# Parametros de entrada
DEVICE=${1:-"default"}
OPTIMIZATION=${2:-"std"}
CPP_FILE=${3:-"tiny_manna.cpp"}

EXE_NAME=$(basename "$CPP_FILE" .cpp)
RESULTS_FILE="$RESULTS_DIR/${DEVICE}_${OPTIMIZATION}.csv"

# Detectar número máximo de hilos
if command -v nproc &> /dev/null; then
  MAX_THREADS=$(nproc)
else
  MAX_THREADS=$(getconf _NPROCESSORS_ONLN)
fi

# Valores de N y numero de ejecuciones
declare -a N_RUNS_ARRAY=(
  "4096:4"
  "8192:4"
  "16384:4"
  "32768:4"
  "65536:4"
)

# Compiladores y flags
declare -a COMPILER_FLAGS_LIST=(
  "g++:-O3 -march=native -funroll-loops -fopenmp"
  "clang++:-O2 -march=native -mavx2 -fopenmp"
  "clang++:-O3 -march=native -funroll-loops -ffast-math -fopenmp"
  "icpx:-xHost -O3 -mavx2 -qopenmp"
)

# Encabezado CSV si no existe
if [ ! -f "$RESULTS_FILE" ]; then
  echo "N,Compiler,Flags,Grains/µs,Total_Time(s),Memory_KB,Threads" \
       > "$RESULTS_FILE"
fi

for N_RUNS in "${N_RUNS_ARRAY[@]}"; do
  N=${N_RUNS%%:*}
  RUNS=${N_RUNS##*:}

  echo -e "\n--- N = $N ($RUNS ejecuciones) ---"

  # Actualizar N en params.h
  if [ -f params.h ]; then
    sed -i "s/^#define N .*/#define N $N/" params.h
  fi

  for ((THREADS=1; THREADS <= MAX_THREADS; THREADS *= 2)); do
    echo " * Probando con THREADS=$THREADS"

    # Actualizar THREADS en params.h
    if [ -f params.h ]; then
      sed -i "s/^#define THREADS .*/#define THREADS $THREADS/" params.h
    fi

    for comp_pair in "${COMPILER_FLAGS_LIST[@]}"; do
      CXX=${comp_pair%%:*}
      FLAGS=${comp_pair#*:}

      if ! command -v "$CXX" &> /dev/null; then
        echo "  - $CXX no disponible, omitiendo..."
        continue
      fi

      # Evitar duplicados
      if grep -q "^${N},${CXX},\"${FLAGS}\",.*,$THREADS\$" "$RESULTS_FILE"; then
        echo "  - Ya registrado para N=$N, $CXX \"$FLAGS\", THREADS=$THREADS. Omitiendo..."
        continue
      fi

      BEST_GRAINS=0
      BEST_TIME=1e9
      BEST_MEM=999999999

      for ((j=1; j<=RUNS; j++)); do
        echo -n "    Ejecutando $j/$RUNS… "

        $CXX $FLAGS "$CPP_FILE" -o "$EXE_NAME" || {
          echo "Falló compilación con $CXX"; break;
        }

        OUTPUT=$(OMP_NUM_THREADS="$THREADS" ./"$EXE_NAME" 2>&1)
        TIME_S=$(echo "$OUTPUT" | awk '/Tiempo de procesamiento/ {print $5}')
        GRAINS=$(echo "$OUTPUT" | awk '/Granos\/us/ {print $2}')
        MEM=$(ps -o rss= -p $$ | tr -d ' ')

        echo "Grains/µs=$GRAINS | Time=${TIME_S}s | Mem=${MEM}KB"

        [[ $GRAINS > $BEST_GRAINS ]] && BEST_GRAINS=$GRAINS
        awk -v t="$TIME_S" -v bt="$BEST_TIME" 'BEGIN{exit !(t<bt)}' && BEST_TIME=$TIME_S
        (( MEM < BEST_MEM )) && BEST_MEM=$MEM
      done

      # Guardar resultados con THREADS al final
      echo "${N},${CXX},\"${FLAGS}\",${BEST_GRAINS},${BEST_TIME},${BEST_MEM},${THREADS}" \
        >> "$RESULTS_FILE"

      echo "    >> Mejor: Grains/µs=$BEST_GRAINS | Time=${BEST_TIME}s | Mem=${BEST_MEM}KB"
    done
  done
done

echo -e "\n=== Benchmarks completados: ==="
column -s, -t "$RESULTS_FILE"

#!/bin/bash

# Directorio para resultados
RESULTS_DIR="benchmarks"
mkdir -p "$RESULTS_DIR"

# Parámetros de entrada
DEVICE=${1:-"default"}        # Identificador de dispositivo
OPTIMIZATION=${2:-"std"}       # Tipo de optimización
CPP_FILE=${3:-"tiny_manna.cpp"} # Archivo .cpp a compilar

# Nombre del ejecutable (sin extensión)
EXE_NAME=$(basename "$CPP_FILE" .cpp)

# Archivo de resultados
RESULTS_FILE="$RESULTS_DIR/${DEVICE}_${OPTIMIZATION}.csv"

# Valores de N y número de ejecuciones
declare -a N_RUNS_ARRAY=(
  "4096:4"
  "8192:4"
  "16384:2"
  "32768:2"
  "65536:1"
  "131072:1"
  "262144:1"
)

# Compiladores y flags
declare -a COMPILER_FLAGS_LIST=(
  "g++:-O3 -march=native -funroll-loops" # -fopenmp"
  "clang++:-O2 -march=native -mavx2" # -fopenmp"
  "clang++:-O3 -march=native -funroll-loops -ffast-math" # -fopenmp"
  "icpx:-xHost -O3 -mavx2" # -qopenmp"
)

# Si el archivo no existe, crear encabezado
if [ ! -f "$RESULTS_FILE" ]; then
  echo "N,Compiler,Flags,Grains_per_us,Total_Time_s,Memory_KB" > "$RESULTS_FILE"
fi

# Bucle principal: para cada N y cada compilador/flags
for N_RUNS in "${N_RUNS_ARRAY[@]}"; do
  N=${N_RUNS%%:*}
  RUNS=${N_RUNS##*:}

  echo -e "\n--- N = $N ($RUNS ejecuciones) ---"

  # Actualizar params.h si existe
  if [ -f params.h ]; then
    sed -i "s/^#define N .*/#define N $N/" params.h
  fi

  for comp_pair in "${COMPILER_FLAGS_LIST[@]}"; do
    CXX=${comp_pair%%:*}
    FLAGS=${comp_pair#*:}

    # Comprobamos compilador
    if ! command -v "$CXX" &> /dev/null; then
      echo "Compilador $CXX no disponible, omitiendo..."
      continue
    fi

    # Verificar si ya existe en CSV
    if grep -q "^${N},${CXX},\"${FLAGS}\"," "$RESULTS_FILE"; then
      echo "Ya hay datos para N=$N, $CXX \"$FLAGS\". Omitiendo..."
      continue
    fi

    # Inicializar mejores métricas
    BEST_GRAINS=0
    BEST_TIME=1e9
    BEST_MEM=999999999

    for ((j=1; j<=RUNS; j++)); do
      echo -n "  Ejecutando $j/$RUNS... "

      # Compilar
      $CXX $FLAGS "$CPP_FILE" -o "$EXE_NAME" || {
        echo "Falló compilación con $CXX"; break;
      }

      # Ejecutar y capturar métricas
      OUTPUT=$(./"$EXE_NAME" 2>&1)
      TIME_S=$(echo "$OUTPUT" | awk '/Tiempo de procesamiento/ {print $5}')
      GRAINS=$(echo "$OUTPUT" | awk '/Granos\/us/ {print $2}')
      MEM=$(ps -o rss= -p $$ | tr -d ' ')

      echo " Grains/µs=$GRAINS | Time=${TIME_S}s | Mem=${MEM}KB"

      # Actualizar resultados óptimos
      awk -v g="$GRAINS" -v bg="$BEST_GRAINS" 'BEGIN{exit !(g>bg)}' && BEST_GRAINS=$GRAINS
      awk -v t="$TIME_S" -v bt="$BEST_TIME" 'BEGIN{exit !(t<bt)}' && BEST_TIME=$TIME_S
      if [ "$MEM" -lt "$BEST_MEM" ]; then
        BEST_MEM=$MEM
      fi
    done

    # Guardar resultado final en CSV
    echo "${N},${CXX},\"${FLAGS}\",${BEST_GRAINS},${BEST_TIME},${BEST_MEM}" \
      >> "$RESULTS_FILE"

    echo ">> Mejor para N=$N, $CXX \"$FLAGS\":"
    echo "   Grains/µs: $BEST_GRAINS | Time: ${BEST_TIME}s | Mem: ${BEST_MEM}KB"
  done
done

# Mostrar resultados finales
echo -e "\n=== Benchmarks completados: ==="
cat "$RESULTS_FILE"
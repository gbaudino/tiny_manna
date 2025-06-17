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

# Detectar si es version 6 o 7 (busca v6 o v7 en el nombre del archivo)
HAS_THREADS=false
if [[ "$CPP_FILE" =~ v[67] ]]; then
    HAS_THREADS=true
    echo "Detectada versión 6 o 7 - usando THREADS"
else
    echo "Versión anterior a v6 - sin soporte para THREADS"
fi

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

# Compiladores y flags según la versión
if [ "$HAS_THREADS" = true ]; then
    # Para versiones v6/v7 con soporte de threads (necesitan OpenMP)
    declare -a COMPILER_FLAGS_LIST=(
        "g++:-O3 -march=native -funroll-loops -fopenmp"
        "clang++:-O2 -march=native -mavx2 -fopenmp -lomp"
        "clang++:-O3 -march=native -funroll-loops -ffast-math -fopenmp -lomp"
        "icpx:-xHost -O3 -mavx2 -qopenmp"
    )
    echo "Usando flags con OpenMP para versión con threads"
else
    # Para versiones anteriores sin threads (sin OpenMP)
    declare -a COMPILER_FLAGS_LIST=(
        "g++:-O3 -march=native -funroll-loops"
        "clang++:-O2 -march=native -mavx2"
        "clang++:-O3 -march=native -funroll-loops -ffast-math"
        "icpx:-xHost -O3 -mavx2"
    )
    echo "Usando flags sin OpenMP para versión sin threads"
fi

# Encabezado CSV si no existe
if [ ! -f "$RESULTS_FILE" ]; then
    if [ "$HAS_THREADS" = true ]; then
        echo "N,Compiler,Flags,Grains/µs,Total_Time(s),Memory_KB,Threads" > "$RESULTS_FILE"
    else
        echo "N,Compiler,Flags,Grains/µs,Total_Time(s),Memory_KB" > "$RESULTS_FILE"
    fi
fi

for N_RUNS in "${N_RUNS_ARRAY[@]}"; do
    N=${N_RUNS%%:*}
    RUNS=${N_RUNS##*:}
    echo -e "\n--- N = $N ($RUNS ejecuciones) ---"
    
    # Actualizar N en params.h
    if [ -f params.h ]; then
        # Buscar si existe el patrón nuevo (constexpr size_t N = ...)
        if grep -q "constexpr size_t N =" params.h; then
            sed -i "s/^constexpr size_t N = [0-9]\+;/constexpr size_t N = $N;/" params.h
        # Buscar si existe el patrón viejo (#define N ...)
        elif grep -q "^#define N " params.h; then
            sed -i "s/^#define N .*/#define N $N/" params.h
        else
            echo "   Advertencia: No se encontró definición de N en params.h"
        fi
    fi
    
    # Configurar bucle de threads según la versión
    if [ "$HAS_THREADS" = true ]; then
        # Para v6/v7: probar diferentes números de threads
        THREAD_VALUES=($(seq 1 2 $MAX_THREADS | head -10))  # 1, 3, 5, 7, ... hasta MAX_THREADS
        THREAD_VALUES+=($(for i in 1 2 4 8 16; do [[ $i -le $MAX_THREADS ]] && echo $i; done))
        # Eliminar duplicados y ordenar
        THREAD_VALUES=($(printf '%s\n' "${THREAD_VALUES[@]}" | sort -nu))
    else
        # Para versiones anteriores: solo una "iteración" sin threads
        THREAD_VALUES=(1)
    fi
    
    for THREADS in "${THREAD_VALUES[@]}"; do
        if [ "$HAS_THREADS" = true ]; then
            echo " * Probando con THREADS=$THREADS"
            # Actualizar THREADS en params.h
            if [ -f params.h ]; then
                # Buscar si existe el patrón nuevo (constexpr size_t THREADS = ...)
                if grep -q "constexpr.*THREADS =" params.h; then
                    sed -i "s/^constexpr.*THREADS = [0-9]\+;/constexpr size_t THREADS = $THREADS;/" params.h
                # Buscar si existe el patrón viejo (#define THREADS ...)
                elif grep -q "^#define THREADS " params.h; then
                    sed -i "s/^#define THREADS .*/#define THREADS $THREADS/" params.h
                else
                    echo "     Advertencia: No se encontró definición de THREADS en params.h"
                fi
            fi
        else
            echo " * Ejecutando sin configuración de THREADS"
        fi
        
        for comp_pair in "${COMPILER_FLAGS_LIST[@]}"; do
            CXX=${comp_pair%%:*}
            FLAGS=${comp_pair#*:}
            
            if ! command -v "$CXX" &> /dev/null; then
                echo "   - $CXX no disponible, omitiendo..."
                continue
            fi
            
            # Evitar duplicados - ajustar según si tiene threads o no
            if [ "$HAS_THREADS" = true ]; then
                GREP_PATTERN="^${N},${CXX},\"${FLAGS}\",.*,$THREADS\$"
            else
                GREP_PATTERN="^${N},${CXX},\"${FLAGS}\",.* [0-9]\+\$"
            fi
            
            if grep -q "$GREP_PATTERN" "$RESULTS_FILE"; then
                if [ "$HAS_THREADS" = true ]; then
                    echo "   - Ya registrado para N=$N, $CXX \"$FLAGS\", THREADS=$THREADS. Omitiendo..."
                else
                    echo "   - Ya registrado para N=$N, $CXX \"$FLAGS\". Omitiendo..."
                fi
                continue
            fi
            
            BEST_GRAINS=0
            BEST_TIME=1e9
            BEST_MEM=999999999
            
            for ((j=1; j<=RUNS; j++)); do
                echo -n "   Ejecutando $j/$RUNS… "
                
                echo "   Compilando: $CXX $FLAGS $CPP_FILE -o $EXE_NAME"
                COMPILE_OUTPUT=$($CXX $FLAGS "$CPP_FILE" -o "$EXE_NAME" 2>&1)
                if [ $? -ne 0 ]; then
                    echo "   ERROR: Falló compilación con $CXX"
                    echo "   Comando: $CXX $FLAGS $CPP_FILE -o $EXE_NAME"
                    echo "   Salida del error:"
                    echo "$COMPILE_OUTPUT" | sed 's/^/     /'
                    break
                fi
                
                # Ejecutar con o sin OMP_NUM_THREADS según la versión
                if [ "$HAS_THREADS" = true ]; then
                    OUTPUT=$(OMP_NUM_THREADS="$THREADS" ./"$EXE_NAME" 2>&1)
                else
                    OUTPUT=$(./"$EXE_NAME" 2>&1)
                fi
                
                TIME_S=$(echo "$OUTPUT" | awk '/Tiempo de procesamiento/ {print $5}')
                GRAINS=$(echo "$OUTPUT" | awk '/Granos\/us/ {print $2}')
                MEM=$(ps -o rss= -p $$ | tr -d ' ')
                
                echo "Grains/µs=$GRAINS | Time=${TIME_S}s | Mem=${MEM}KB"
                
                [[ $GRAINS > $BEST_GRAINS ]] && BEST_GRAINS=$GRAINS
                awk -v t="$TIME_S" -v bt="$BEST_TIME" 'BEGIN{exit !(t<bt)}' && BEST_TIME=$TIME_S
                [[ $MEM -lt $BEST_MEM ]] && BEST_MEM=$MEM
            done
            
            # Escribir resultado al CSV según la versión
            if [ "$HAS_THREADS" = true ]; then
                echo "$N,$CXX,\"$FLAGS\",$BEST_GRAINS,$BEST_TIME,$BEST_MEM,$THREADS" >> "$RESULTS_FILE"
                echo "   >> Mejor: Grains/µs=$BEST_GRAINS | Time=${BEST_TIME}s | Mem=${BEST_MEM}KB | Threads=$THREADS"
            else
                echo "$N,$CXX,\"$FLAGS\",$BEST_GRAINS,$BEST_TIME,$BEST_MEM" >> "$RESULTS_FILE"
                echo "   >> Mejor: Grains/µs=$BEST_GRAINS | Time=${BEST_TIME}s | Mem=${BEST_MEM}KB"
            fi
        done
    done
done

echo -e "\n=== Benchmarks completados ==="
column -s, -t "$RESULTS_FILE"
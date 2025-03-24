#!/bin/bash

# Número de ejecuciones por prueba
RUNS=${1:-5}

# Directorio para resultados
RESULTS_DIR="benchmarks"
mkdir -p $RESULTS_DIR

# Archivo de resultados
RESULTS_FILE="$RESULTS_DIR/$(date +%m%d_%H%M).csv"

# Encabezado del archivo CSV
echo "Compilador,Flags,Mejor_Metrica" > $RESULTS_FILE

# Lista de compiladores a probar
COMPILERS=(
    "g++"
    "clang++"
    #"icpx"
    #"nvc++"
)

# Flags de compilación a probar
FLAGS_LIST=(
    "-O3 -march=native"
    "-O3 -march=native -funroll-loops"
    "-O3 -march=native -funroll-loops -ffast-math"
    "-O3 -march=native -funroll-loops -ffast-math -ftree-vectorize"
    "-O3 -march=native -funroll-loops -ffast-math -ftree-vectorize -flto"
)

for CXX in "${COMPILERS[@]}"; do
    # Verificar si el compilador está instalado
    if ! command -v $CXX &> /dev/null; then
        echo "Compilador $CXX no encontrado, omitiendo..."
        continue
    fi
    
    for FLAGS in "${FLAGS_LIST[@]}"; do
        echo "Probando con $CXX y flags: $FLAGS"
        
        # Limpiar y compilar
        make clean > /dev/null 2>&1
        OPTFLAGS="$FLAGS" CXX=$CXX make > /dev/null 2>&1

        if [ ! -f tiny_manna ]; then
            echo "Error compilando con $CXX y flags $FLAGS, omitiendo..."
            continue
        fi
        
        # Ejecutar el benchmark
        echo "Ejecutando $RUNS veces..."
        BEST_GRANOS=0
        
        for ((i=1; i<=$RUNS; i++)); do
            echo -n "Ejecución $i/$RUNS... "
            
            # Ejecutar y obtener granos/µs directamente
            OUTPUT=$(./tiny_manna)
            
            # Extraer el valor de granos/µs (asumiendo que es la última línea o valor numérico)
            GRANOS=$(echo "$OUTPUT" | grep -o '[0-9.]*$')
            
            echo "Granos/µs: $GRANOS"
            
            # Verificar si es el mejor resultado
            if (( $(echo "$GRANOS > $BEST_GRANOS" | bc -l) )); then
                BEST_GRANOS=$GRANOS
            fi
        done
        
        # Guardar resultados
        echo "$CXX,$FLAGS,$BEST_GRANOS" >> $RESULTS_FILE
        
        # Mostrar resultados parciales
        echo "Mejor resultado para $CXX con flags $FLAGS:"
        echo "  Granos/µs: $BEST_GRANOS"
        echo "-----------------------------------------"
    done
done
# Ultimo clean
make clean > /dev/null 2>&1
# Ordenar resultados por granos/µs (mejor a peor)
echo -e "\nResultados ordenados por rendimiento (mejor a peor):"
# TODO: fix sorting (this consider the csv header as a row)
sort -t, -k3 -nr $RESULTS_FILE

echo -e "\nTodos los resultados guardados en $RESULTS_FILE"
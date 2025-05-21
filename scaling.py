#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt

def main(csv_path):
    # 1) Carga de datos
    df = pd.read_csv(csv_path)

    # 2) Para cada (Threads, N) quedarnos con el mejor Grains/µs
    best = df.groupby(['Threads', 'N'])['Grains/µs'] \
             .max().reset_index()

    # 3) Pivot: Threads en filas, N en columnas
    pivot = best.pivot(index='Threads', columns='N', values='Grains/µs')


    # 5) Graficar curvas: eje X = Threads, eje Y = Grains/µs
    plt.figure(figsize=(10,6))
    for N in sorted(pivot.columns):
        plt.plot(pivot.index, pivot[N],
                 marker='o', label=f'N={N}')
    plt.xlabel('Threads')
    plt.ylabel('Grains/µs')
    plt.title('Escalabilidad de Threads')
    plt.xscale('log', base=2)
    plt.xticks([1, 2, 4, 8], [1, 2, 4, 8])
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.legend(title='Tamaño N', reverse=True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python plot_scaling.py <ruta_al_csv>")
        sys.exit(1)
    main(sys.argv[1])

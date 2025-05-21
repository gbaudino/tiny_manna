import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from tabulate import tabulate

# Definir colores UNCC
colors = {
    'UnccGreen': '#005035',
    'UnccLightGreen': '#C3D7A4',
    'UnccGold': '#A49665',
    'UnccOrange': '#F3901D',
    'UnccLightYellow': '#899064',
    'UnccBlue': '#007377',
    'UnccPink': '#DE3A6E',
    'UnccPurple': '#A65E9E',
    'White': '#FFFFFF',
    'LightGray': '#F1E6B2'
}

# Estilo de gráficos personalizado
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

def load_best_by_N(filename, is_base=False):
    df = pd.read_csv(filename)
    
    if is_base:
        if 'Flags' in df.columns and 'Compiler' in df.columns:
            df = df[
                (df['Compiler'].str.lower().str.strip() == 'g++') &
                (
                    df['Flags'].isnull() |
                    (df['Flags'].str.strip() == '') |
                    (df['Flags'].str.lower().str.strip() == 'none')
                )
            ]
        else:
            raise ValueError(f"El archivo {filename} necesita columnas 'flags' y 'compiler'.")
        
        if df.empty:
            raise ValueError(f"No se encontraron filas base (g++ sin flags) en {filename}")
        
        base = df.groupby('N')['Grains/µs'].mean().rename(filename)
        return base
    else:
        best = df.groupby('N')['Grains/µs'].max().rename(filename)
        return best

def geometric_mean(lst):
    return np.exp(np.mean(np.log(lst)))

def main(files):
    best_by_file = [load_best_by_N(f, is_base=(i == 0)) for i, f in enumerate(files)]
    all_N = sorted(set().union(*[b.index for b in best_by_file]))
    data = pd.concat(best_by_file, axis=1).loc[all_N]
    data.columns = [f'v{i}' for i in range(len(files))]

    speedup_table = data.copy()
    for i in range(1, len(files)):
        speedup_table[f'v{i-1} → v{i}'] = data.iloc[:, i] / data.iloc[:, i - 1]

    print(tabulate(speedup_table.reset_index(), headers='keys', tablefmt='github', floatfmt=".2f"))

    for i in range(1, len(files)):
        speedups = data.iloc[:, i] / data.iloc[:, i - 1]
        print(f"\nSpeedup general de v{i-1} a v{i}: {geometric_mean(speedups):.2f}×")

    plt.figure()
    color_list = list(colors.values())
    for i, col in enumerate(data.columns[:len(files)]):
        plt.plot(data.index, data[col], marker='o', label=col, color=color_list[i % len(color_list)])

    plt.xlabel('N')
    plt.ylabel('Mejor granos/us')
    plt.title('Comparación de rendimiento por versión')
    plt.legend()
    plt.xscale('log', base=2)
    plt.grid(True)
    # Set custom background color
    bg_color = '#FAFAFA'
    plt.gcf().patch.set_facecolor(bg_color)
    plt.gca().set_facecolor(bg_color)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python compare_speedup.py archivo1.csv archivo2.csv ...")
    else:
        main(sys.argv[1:])

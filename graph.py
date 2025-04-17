import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Definir los colores UNCC
colors = {
    'UnccGreen': '#005035',
    'UnccLightGreen': '#C3D7A4',
    'UnccGold': '#A49665',
    'UnccOrange': '#F3901D',
    'UnccLightYellow': '#899064',
    'UnccBlue': '#007377',
    'UnccPink': '#DE3A6E',
    'White': '#FFFFFF',
    'LightGray': '#F1E6B2'
}

# Crear un estilo personalizado para los gráficos
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

# Cargar los datos originales (desktop)
v0_df = pd.read_csv('benchmarks/desktop_v0.csv')
v1_df = pd.read_csv('benchmarks/desktop_v1.csv')
v2_df = pd.read_csv('benchmarks/desktop_v2.csv')

# Cargar los datos nuevos (M3) - solo para comparaciones específicas
m3_v0_df = pd.read_csv('benchmarks/M3_v0.csv')
m3_v1_df = pd.read_csv('benchmarks/M3_v1.csv')
m3_v2_df = pd.read_csv('benchmarks/M3_v2.csv')

# Agregar columnas para identificar las versiones
v0_df['approach'] = 'Version 0'
v1_df['approach'] = 'Version 1'
v2_df['approach'] = 'Version 2'
v0_df['platform'] = 'Desktop'
v1_df['platform'] = 'Desktop'
v2_df['platform'] = 'Desktop'

m3_v0_df['approach'] = 'Version 0'
m3_v1_df['approach'] = 'Version 1'
m3_v2_df['approach'] = 'Version 2'
m3_v0_df['platform'] = 'M3'
m3_v1_df['platform'] = 'M3'
m3_v2_df['platform'] = 'M3'

# Combinar los datos para mejor análisis
desktop_data = pd.concat([v0_df, v1_df, v2_df], ignore_index=True)
m3_data = pd.concat([m3_v0_df, m3_v1_df, m3_v2_df], ignore_index=True)
all_data = pd.concat([desktop_data, m3_data], ignore_index=True)

# Función para encontrar el mejor flag para cada compilador en cada enfoque
def find_best_flags(df, approach):
    best_flags = {}
    approach_data = df[df['approach'] == approach]
    
    for compiler in approach_data['Compiler'].unique():
        compiler_data = approach_data[approach_data['Compiler'] == compiler]
        # Agrupar por N y encontrar el mejor flag para cada N
        for n in compiler_data['N'].unique():
            n_data = compiler_data[compiler_data['N'] == n]
            best_idx = n_data['Grains/µs'].idxmax()
            best_flag = n_data.loc[best_idx, 'Flags']
            if compiler not in best_flags:
                best_flags[compiler] = {}
            best_flags[compiler][n] = best_flag
    return best_flags

# Mapeo de compiladores a colores
compiler_colors = {
    'g++': colors['UnccGreen'],
    'clang++': colors['UnccOrange'],
    'icpx': colors['UnccBlue'],
    'nvc++': colors['UnccPink']
}

# Mapeo de enfoques a marcadores
approach_markers = {
    'Version 0': '^',
    'Version 1': 'o',
    'Version 2': 's'
}

# Mapeo de plataformas a estilos de línea
platform_styles = {
    'Desktop': '-',
    'M3': '--'
}

compilers = ['g++', 'clang++', 'icpx', 'nvc++']
m3_compilers = ['g++', 'clang++']  # M3 solo tiene estos compiladores
approaches = ['Version 0', 'Version 1', 'Version 2']
N_values = sorted(desktop_data['N'].unique())

# Colores para los diferentes enfoques
approach_colors = {
    'Version 0': colors['UnccBlue'],
    'Version 1': colors['UnccGold'],
    'Version 2': colors['UnccOrange']
}

# 1. Gráfico de barras: Base vs Optimizado por compilador
def plot_base_vs_optimized(desktop_data):
    plt.figure(figsize=(14, 8))
    
    # Encontrar los mejores flags para cada versión
    best_flags_v0 = find_best_flags(desktop_data, 'Version 0')
    best_flags_v1 = find_best_flags(desktop_data, 'Version 1')
    best_flags_v2 = find_best_flags(desktop_data, 'Version 2')
    
    # Configurar el ancho de las barras
    bar_width = 0.2
    index = np.arange(len(compilers))
    
    for i, approach in enumerate(approaches):
        base_perf = []
        optimized_perf = []
        
        for compiler in compilers:
            # Obtener datos base (sin flags)
            base_data = desktop_data[(desktop_data['approach'] == approach) & 
                                   (desktop_data['Compiler'] == compiler) & 
                                   (desktop_data['Flags'].isin(["", " ", np.nan]))]
            # Obtener datos optimizados
            if approach == 'Version 0':
                best_flag = best_flags_v0[compiler][max(best_flags_v0[compiler].keys())]
            elif approach == 'Version 1':
                best_flag = best_flags_v1[compiler][max(best_flags_v1[compiler].keys())]
            else:
                best_flag = best_flags_v2[compiler][max(best_flags_v2[compiler].keys())]
                
            optimized_data = desktop_data[(desktop_data['approach'] == approach) & 
                                        (desktop_data['Compiler'] == compiler) & 
                                        (desktop_data['Flags'] == best_flag)]
            
            # Usar el promedio de rendimiento para N=4096 como representativo
            base_avg = base_data[base_data['N'] == 4096]['Grains/µs'].mean()
            optimized_avg = optimized_data[optimized_data['N'] == 4096]['Grains/µs'].mean()
            
            base_perf.append(base_avg)
            optimized_perf.append(optimized_avg)
        
        # Dibujar barras
        plt.bar(index + i*bar_width, base_perf, bar_width, 
                color=approach_colors[approach], alpha=0.6, 
                label=f'{approach} Base', hatch='/')
        plt.bar(index + i*bar_width, optimized_perf, bar_width, 
                color=approach_colors[approach], alpha=0.9, 
                label=f'{approach} Optimized', bottom=base_perf)
    
    plt.xlabel('Compiler')
    plt.ylabel('Performance (Grains/µs)')
    plt.title('Base vs Optimized Performance by Compiler and Approach (N=4096)')
    plt.xticks(index + bar_width, compilers)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('base_vs_optimized.png', dpi=300, bbox_inches='tight')
    plt.show()

# 2. Comparación entre plataformas (Desktop vs M3)
def plot_platform_comparison(all_data):
    plt.figure(figsize=(12, 8))
    
    # Seleccionar los mejores compiladores de cada plataforma
    desktop_best = all_data[(all_data['platform'] == 'Desktop') & 
                          (all_data['Compiler'] == 'icpx')]  # icpx fue generalmente el mejor en Desktop
    m3_best = all_data[(all_data['platform'] == 'M3') & 
                     (all_data['Compiler'] == 'clang++')]  # clang++ fue generalmente el mejor en M3
    
    for approach in approaches:
        # Filtrar por enfoque y obtener los mejores flags
        desktop_approach = desktop_best[desktop_best['approach'] == approach]
        m3_approach = m3_best[m3_best['approach'] == approach]
        
        # Obtener los mejores flags para cada N
        desktop_perf = []
        m3_perf = []
        common_Ns = sorted(set(desktop_approach['N']).intersection(set(m3_approach['N'])))
        
        for n in common_Ns:
            desktop_n = desktop_approach[desktop_approach['N'] == n]
            m3_n = m3_approach[m3_approach['N'] == n]
            
            desktop_best_flag = desktop_n.loc[desktop_n['Grains/µs'].idxmax(), 'Flags']
            m3_best_flag = m3_n.loc[m3_n['Grains/µs'].idxmax(), 'Flags']
            
            desktop_perf.append(desktop_n[desktop_n['Flags'] == desktop_best_flag]['Grains/µs'].values[0])
            m3_perf.append(m3_n[m3_n['Flags'] == m3_best_flag]['Grains/µs'].values[0])
        
        plt.plot(common_Ns, desktop_perf, marker=approach_markers[approach], 
                color=approach_colors[approach], linestyle='-',
                label=f'Desktop (icpx) - {approach}')
        plt.plot(common_Ns, m3_perf, marker=approach_markers[approach], 
                color=approach_colors[approach], linestyle='--',
                label=f'M3 (clang++) - {approach}')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Problem Size (N)')
    plt.ylabel('Performance (Grains/µs)')
    plt.title('Performance Comparison: Desktop vs M3 (Best Compilers)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    plt.savefig('platform_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

# 3. Heatmap de speedup por versión y tamaño de problema
def plot_speedup_heatmap(desktop_data):
    plt.figure(figsize=(12, 8))
    
    # Calcular speedup respecto a Version 0 para cada N y compilador
    speedup_data = []
    compilers = ['g++', 'clang++', 'icpx', 'nvc++']
    Ns = sorted(desktop_data['N'].unique())
    
    for n in Ns:
        row = []
        for compiler in compilers:
            # Obtener mejor performance para cada versión
            v0_perf = desktop_data[(desktop_data['N'] == n) & 
                                 (desktop_data['Compiler'] == compiler) & 
                                 (desktop_data['approach'] == 'Version 0')]['Grains/µs'].max()
            v1_perf = desktop_data[(desktop_data['N'] == n) & 
                                 (desktop_data['Compiler'] == compiler) & 
                                 (desktop_data['approach'] == 'Version 1')]['Grains/µs'].max()
            v2_perf = desktop_data[(desktop_data['N'] == n) & 
                                 (desktop_data['Compiler'] == compiler) & 
                                 (desktop_data['approach'] == 'Version 2')]['Grains/µs'].max()
            
            # Calcular speedup relativo a Version 0
            speedup_v1 = v1_perf / v0_perf if v0_perf > 0 else 1
            speedup_v2 = v2_perf / v0_perf if v0_perf > 0 else 1
            row.extend([speedup_v1, speedup_v2])
        speedup_data.append(row)
    
    # Crear el heatmap
    speedup_array = np.array(speedup_data)
    plt.imshow(speedup_array, cmap='RdYlGn', aspect='auto', vmin=0.5, vmax=2.0)
    
    # Configurar ejes
    plt.xticks(np.arange(len(compilers)*2), 
              [f'{comp}\nV1' for comp in compilers] + [f'{comp}\nV2' for comp in compilers])
    plt.yticks(np.arange(len(Ns)), Ns)
    plt.xlabel('Compiler and Version')
    plt.ylabel('Problem Size (N)')
    plt.title('Speedup Relative to Version 0')
    
    # Añadir barra de color
    cbar = plt.colorbar()
    cbar.set_label('Speedup (x times)')
    
    # Añadir valores en las celdas
    for i in range(len(Ns)):
        for j in range(len(compilers)*2):
            plt.text(j, i, f"{speedup_array[i, j]:.1f}x", 
                    ha="center", va="center", color="black", fontsize=8)
    
    plt.tight_layout()
    plt.savefig('speedup_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

# Generar todos los gráficos
plot_base_vs_optimized(desktop_data)
plot_platform_comparison(all_data)
plot_speedup_heatmap(desktop_data)
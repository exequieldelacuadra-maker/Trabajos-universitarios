# ==========================================================================
# 1. IMPORTAR LIBRERÍAS
# ==========================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skfda  # Equivalente a 'fda' en R
from sklearn.mixture import GaussianMixture # Equivalente a 'mclust' en R
import warnings

# Opcional: Ignorar advertencias numéricas
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==========================================================================
# 2. CARGA Y LIMPIEZA DE DATOS (Equivalente a dplyr/tidyr)
# ==========================================================================

# Cargar el archivo
try:
    df_raw = pd.read_csv("huella_carbono_PowerQuery.xlsx - Hoja1.csv")
except FileNotFoundError:
    print("Error: No se encontró el archivo CSV. Asegúrate de que esté en el directorio correcto.")
    # Salir o manejar el error

# --- Limpieza datos (Equivalente a mutate y filter) ---
df_filtrado = df_raw.copy()
df_filtrado['Country'] = df_filtrado['Country'].str.strip().str.replace(r'\\n', '', regex=False)
df_filtrado = df_filtrado[
    df_filtrado['ISO 3166-1 alpha-3'].notna() &
    (df_filtrado['Country'] != 'Global')
]

# --- Pivotear datos (Equivalente a pivot_wider) ---
df_wide = df_filtrado.pivot(
    index='Year',
    columns='Country',
    values='Total'
)

# --- Reemplazar NA con 0 ---
df_wide = df_wide.fillna(0)

# --- Filtrar años (¡NUEVO FILTRO!) ---
df_wide.index = pd.to_numeric(df_wide.index)
# Filtramos de 1900 a 2024
datos_limpios = df_wide[
    (df_wide.index >= 1900) & (df_wide.index <= 2024)
].copy()

# --- Preparar variables para FDA ---
# tiempo_0 <- datos_limpios$Year
tiempo_0 = datos_limpios.index.to_numpy()
print(f"Clase de dato 'tiempo': {type(tiempo_0)}")

# --- Crear matriz numérica final ---
# Esto reemplaza 'datos_numericos' y 'as.matrix'
variables_matriz = datos_limpios.to_numpy()

print(f"Clase de variables: {type(variables_matriz)}")

# Guardamos los nombres de los países
nombres_paises = datos_limpios.columns.to_numpy()

# --- Asegurarnos que es numérico (float/double) ---
variables_matriz = variables_matriz.astype(float)
print(f"Dimensiones de variables_matriz: {variables_matriz.shape}")

# --- Limpieza de valores no finitos (NA, NaN, Inf) ---
if not np.all(np.isfinite(variables_matriz)):
    print("Encontrados valores no finitos (NaN, Inf)...")
    # Reemplazar no-finitos con 0 (como en el script)
    variables_matriz[~np.isfinite(variables_matriz)] = 0

# (Tu script R también reemplazaba por la media, pero
#  primero reemplazaba por 0. Este paso ya es muy robusto.)
# Opcional: Reemplazar 0s por la media de la columna si se desea.

# ==========================================================================
# 3. VISUALIZACIÓN GENERAL (Equivalente a matplot)
# ==========================================================================

# --- Muestra de 30 países ---
# np.random.seed(123) # Opcional: para reproducibilidad
indices_muestra = np.random.choice(variables_matriz.shape[1], 30, replace=False)

plt.figure(figsize=(10, 6))
plt.plot(tiempo_0, variables_matriz[:, indices_muestra])
plt.xlabel("Año")
plt.ylabel("Emisiones Totales")
plt.title("Curvas de Emisión (Muestra de 30 países aleatorios)")
plt.show()

# --- Muestra de 10 países ---
indices_muestra_1 = np.random.choice(variables_matriz.shape[1], 10, replace=False)
plt.figure(figsize=(10, 6))
plt.plot(tiempo_0, variables_matriz[:, indices_muestra_1])
plt.xlabel("Año")
plt.ylabel("Emisiones Totales")
plt.title("Curvas de Emisión (Muestra de 10 países aleatorios)")
plt.show()

# --- Muestra de 5 países ---
indices_muestra_2 = np.random.choice(variables_matriz.shape[1], 5, replace=False)
plt.figure(figsize=(10, 6))
plt.plot(tiempo_0, variables_matriz[:, indices_muestra_2])
plt.xlabel("Año")
plt.ylabel("Emisiones Totales")
plt.title("Curvas de Emisión (Muestra de 5 países aleatorios)")
plt.show()

# --- Muestra de 1 país (con leyenda) ---
indices_muestra_3 = np.random.choice(variables_matriz.shape[1], 1, replace=False)
# Obtener el nombre
nombre_pais = nombres_paises[indices_muestra_3[0]]

plt.figure(figsize=(10, 6))
plt.plot(tiempo_0, variables_matriz[:, indices_muestra_3], label=nombre_pais, color='blue', linewidth=2)
plt.xlabel("Año")
plt.ylabel("Emisiones Totales")
plt.title(f"Curva de Emisión de: {nombre_pais}")
plt.legend() # <-- Añade la leyenda
plt.show()


# ==========================================================================
# 4. ANÁLISIS FDA (Equivalente a fda con scikit-fda)
# ==========================================================================

# --- Actualización variables ---
tiempo = tiempo_0
# ¡¡APLICANDO LOGARITMO!!
variables = np.log(variables_matriz + 1) 
# (nombres_paises ya está definido)


# --- Crear FDataGrid (el objeto de datos funcionales en skfda) ---
# NOTA: skfda espera los datos como (n_muestras, n_puntos_de_tiempo)
# Por lo tanto, debemos TRANSPONER nuestra matriz 'variables'.
fd_grid = skfda.FDataGrid(
    data_matrix=variables.T, # Transponer!
    grid_points=tiempo
)

# --- Creación de las bases (nbasis=10) ---
rango_tiempo = (tiempo.min(), tiempo.max())
nbasis_sugerido = 10  # <-- CAMBIO DE TU SCRIPT

base_bspline = skfda.representation.basis.BSpline(
    domain_range=rango_tiempo, 
    n_basis=nbasis_sugerido
)
base_fourier = skfda.representation.basis.Fourier(
    domain_range=rango_tiempo,
    n_basis=nbasis_sugerido # n_basis=10 es 1 (const) + 4 pares sin/cos + 1 sin extra
)

# --- Ajustar coeficientes (Equivalente a Data2fd) ---
repre_bspline = fd_grid.to_basis(base_bspline)
repre_fourier = fd_grid.to_basis(base_fourier)

# --- Gráficas (Equivalente a plot.fd) ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

repre_bspline.plot(axes=axes[0])
axes[0].set_title("Gráfica DF - B-Spline (Log)")
axes[0].set_xlabel("Año")
axes[0].set_ylabel("Log(Emisiones)")

repre_fourier.plot(axes=axes[1])
axes[1].set_title("Gráfica DF - Fourier (Log)")
axes[1].set_xlabel("Año")
axes[1].set_ylabel("Log(Emisiones)")

plt.tight_layout()
plt.show()

# --- Evaluación del Error (MSE) ---
# Usamos 'variables' (la matriz logarítmica) como 'y'
y_pred_bspline = repre_bspline.evaluate(tiempo).T 
y_pred_fourier = repre_fourier.evaluate(tiempo).T

mse_bspline = np.mean((variables - y_pred_bspline)**2)
mse_fourier = np.mean((variables - y_pred_fourier)**2)

print(f"MSE B-Spline: {mse_bspline:.4f}")
print(f"MSE Fourier: {mse_fourier:.4f}")

# --- Media y Varianza Funcional (B-SPLINE) ---
mean_fd = repre_bspline.mean()
cov_fd = repre_bspline.cov()

# --- Gráfico de curvas, media y varianza ---
grid_time = np.linspace(rango_tiempo[0], rango_tiempo[1], 200)

fd_values_bspline = repre_bspline.evaluate(grid_time).T # (n_puntos, n_curvas)
mean_values_bspline = mean_fd.evaluate(grid_time).flatten() # vector 1D
var_values_bspline = np.diag(cov_fd.evaluate(grid_time, grid_time))
var_values_bspline[var_values_bspline < 0] = 0 # Corregir errores numéricos
sd_values_bspline = np.sqrt(var_values_bspline)

plt.figure(figsize=(10, 6))
# matplot(..., col="lightgray")
plt.plot(grid_time, fd_values_bspline, color='lightgray', alpha=0.5)
# lines(..., col="blue")
plt.plot(grid_time, mean_values_bspline, color='blue', linewidth=3, label='Media')
# lines(..., col="red")
plt.plot(grid_time, mean_values_bspline + sd_values_bspline, color='red', linestyle='--', linewidth=2, label='Media ± 1 SD')
plt.plot(grid_time, mean_values_bspline - sd_values_bspline, color='red', linestyle='--', linewidth=2)

plt.legend(loc="upper left")
plt.xlabel("Año")
plt.ylabel("Log(Emisiones)") # Etiqueta correcta
plt.title("B-Spline: Curvas de Países, Media y Varianza (Log)")
plt.show()

# ==========================================================================
# 5. FPCA (Análisis de Componentes Principales Funcionales)
# ==========================================================================

# Usamos B-Spline
nharm_max = nbasis_sugerido - 1 # 10 - 1 = 9
fpca_bspline = skfda.preprocessing.dim_reduction.FPCA(
    n_components=nharm_max
)
fpca_bspline.fit(repre_bspline)

# Varianza acumulada
cumvar_bspline = np.cumsum(fpca_bspline.explained_variance_ratio_)
try:
    k95_bspline = np.where(cumvar_bspline >= 0.95)[0][0] + 1 # +1 por índice 0
    print(f"B-Spline necesita {k95_bspline} componentes para ≥95% de varianza")
except IndexError:
    print(f"No se alcanzó el 95% de varianza con {nharm_max} componentes.")
    k95_bspline = nharm_max

# Graficar autofunciones
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fpca_bspline.components_[0].plot(axes=axes[0], color='blue')
axes[0].set_title("B-Spline - Autofunción CP1") # Título añadido
axes[0].set_xlabel("Año")

fpca_bspline.components_[1].plot(axes=axes[1], color='red')
axes[1].set_title("B-Spline - Autofunción CP2") # Título añadido
axes[1].set_xlabel("Año")
plt.show()

# ==========================================================================
# 6. CLUSTERING FUNCIONAL (Equivalente a mclust)
# ==========================================================================

# Obtener Scores
scores = fpca_bspline.transform(repre_bspline)
scores_bspline = scores[:, 0:2] # Tomamos los 2 primeros

# Ajustar modelo (GaussianMixture)
n_components_range = range(2, 7) # R's 2:6
models = [
    GaussianMixture(n, covariance_type='full', random_state=123).fit(scores_bspline)
    for n in n_components_range
]
bics = [m.bic(scores_bspline) for m in models]
best_model_index = np.argmin(bics)
mclust_result = models[best_model_index]
n_clusters_optimo = n_components_range[best_model_index]
cluster_asignado = mclust_result.predict(scores_bspline)

print(f"El número óptimo de clusters (por BIC) es: {n_clusters_optimo}")

# Graficar los clusters
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    scores_bspline[:, 0], 
    scores_bspline[:, 1], 
    c=cluster_asignado, 
    cmap='rainbow', 
    alpha=0.7
)
plt.xlabel("Score Componente Principal 1")
plt.ylabel("Score Componente Principal 2")
plt.title("Clusters de Países (Gaussian Mixture)")
plt.colorbar(scatter, label="ID de Clúster")
plt.show()

# ==========================================================================
# 7. VISUALIZACIÓN FINAL POR CLUSTER
# ==========================================================================

# Crear un dataframe para ver qué país quedó en qué cluster
df_clusters_paises = pd.DataFrame({
    'Pais': nombres_paises,
    'Cluster': cluster_asignado,
    'Score_CP1': scores_bspline[:, 0],
    'Score_CP2': scores_bspline[:, 1]
})

print("\n--- Países por cluster ---")
for k in range(n_clusters_optimo):
    paises_en_cluster = df_clusters_paises[df_clusters_paises['Cluster'] == k]['Pais'].values
    print(f"\nCluster {k} ({len(paises_en_cluster)} países):")
    print(", ".join(paises_en_cluster))

# --- Graficar las curvas originales coloreadas por el cluster ---
plt.figure(figsize=(12, 8))
colors = plt.cm.rainbow(np.linspace(0, 1, n_clusters_optimo))

for k in range(n_clusters_optimo):
    indices_cluster = np.where(cluster_asignado == k)[0]
    plt.plot(grid_time, fd_values_bspline[:, indices_cluster], color=colors[k], alpha=0.2)

# --- Agregar la media funcional de cada clúster ---
for k in range(n_clusters_optimo):
    indices_cluster = np.where(cluster_asignado == k)[0]
    mean_cluster = repre_bspline[indices_cluster].mean()
    mean_values = mean_cluster.evaluate(grid_time).flatten()
    plt.plot(grid_time, mean_values, color=colors[k], linewidth=4, linestyle='--', label=f'Media Clúster {k}')

plt.xlabel("Año")
plt.ylabel("Log(Emisiones)")
plt.title("Curvas de Emisión Coloreadas por Clúster (Log)")
plt.legend(loc="upper left")
plt.show()

# --- Gráficas separadas por clúster ---
print("\n--- Generando gráficos separados por clúster ---")
for k in range(n_clusters_optimo):
    fig, ax = plt.subplots(figsize=(10, 6))
    indices_cluster = np.where(cluster_asignado == k)[0]
    
    if len(indices_cluster) == 0:
        continue # Omitir clusters vacíos si los hubiera

    # 'curvas_cluster'
    curvas_cluster = fd_values_bspline[:, indices_cluster]
    # 'mean_cluster'
    mean_cluster = repre_bspline[indices_cluster].mean()
    mean_values = mean_cluster.evaluate(grid_time).flatten()

    ax.plot(grid_time, curvas_cluster, color=colors[k], alpha=0.4)
    ax.plot(grid_time, mean_values, color=colors[k], linewidth=4, label=f'Media funcional')
    
    ax.set_xlabel("Año")
    ax.set_ylabel("Log(Emisiones)")
    ax.set_title(f"Curvas de Emisión - Clúster {k} ({len(indices_cluster)} países)")
    ax.legend(loc="upper left")
    plt.show()
    
    # Mostrar nombres de países en consola
    print(f"\n====================\nCluster {k} contiene los países:\n")
    print(", ".join(df_clusters_paises[df_clusters_paises['Cluster'] == k]['Pais'].values))
    print("====================\n")

# --- Gráfico separado para el Clúster 0 (Equivalente a k=1 en R) ---
print("\n--- Gráfico solo para el primer clúster (k=0) ---")
k = 0 # Usamos k=0 para el primer clúster en Python
fig, ax = plt.subplots(figsize=(10, 6))
indices_cluster = np.where(cluster_asignado == k)[0]

if len(indices_cluster) > 0:
    curvas_cluster = fd_values_bspline[:, indices_cluster]
    mean_cluster = repre_bspline[indices_cluster].mean()
    mean_values = mean_cluster.evaluate(grid_time).flatten()

    ax.plot(grid_time, curvas_cluster, color=colors[k], alpha=0.4)
    ax.plot(grid_time, mean_values, color=colors[k], linewidth=4, label='Media funcional')
    
    ax.set_xlabel("Año")
    ax.set_ylabel("Log(Emisiones)")
    ax.set_title(f"Curvas de Emisión - Clúster {k} ({len(indices_cluster)} países)")
    ax.legend(loc="upper left")
    plt.show()
else:
    print(f"El Clúster {k} está vacío.")

print("\n--- Traducción a Python completada ---")
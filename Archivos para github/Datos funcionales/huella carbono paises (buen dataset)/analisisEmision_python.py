# ==========================================================================
# 1. IMPORTAR LIBRERÍAS
# ==========================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skfda  # Equivalente a 'fda' en R
from sklearn.mixture import GaussianMixture # Equivalente a 'mclust' en R
import warnings

# Opcional: Ignorar advertencias numéricas comunes de skfda
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==========================================================================
# 2. CARGA Y LIMPIEZA DE DATOS (Equivalente a dplyr/tidyr)
# ==========================================================================

# Cargar el archivo
try:
    # Asumimos que el CSV está en el mismo directorio
    df_raw = pd.read_excel("huella_carbono_PowerQuery.xlsx")
except FileNotFoundError:
    print("Error: No se encontró el archivo CSV. Asegúrate de que esté en el directorio correcto.")
    # Salir o manejar el error
    
# print(df_raw.head()) # Equivalente a head()

# --- Limpieza datos (Equivalente a mutate y filter) ---
df_filtrado = df_raw.copy()
# Equivalente a str_trim y str_replace_all
df_filtrado['Country'] = df_filtrado['Country'].str.strip().str.replace(r'\\n', '', regex=False)

# Equivalente a filter(!is.na(...) & Country != "Global")
df_filtrado = df_filtrado[
    df_filtrado['ISO 3166-1 alpha-3'].notna() &
    (df_filtrado['Country'] != 'Global')
]

# --- Pivotear datos (Equivalente a pivot_wider) ---
# En pandas, 'pivot' hace esto. El índice se convierte en 'id_cols' (Year)
df_wide = df_filtrado.pivot(
    index='Year',
    columns='Country',
    values='Total'
)

# --- Reemplazar NA con 0 ---
df_wide = df_wide.fillna(0)

# --- Filtrar fila de "basura" (el año 535 que encontraste) ---
# Asegurarnos de que el índice (Year) sea numérico
df_wide.index = pd.to_numeric(df_wide.index)
df_limpios = df_wide[
    (df_wide.index >= 1750) & (df_wide.index <= 2024)
].copy()

# --- Preparar variables para FDA ---
# Equivalente a: tiempo_0 <- datos_limpios$Year
tiempo_0 = df_limpios.index.to_numpy()

# Equivalente a todos los pasos de 'as.matrix', 'mutate(across...)'
# .to_numpy() crea la matriz numérica directamente.
variables_matriz = df_limpios.to_numpy()

# Limpieza de valores no finitos (NaN, Inf)
if not np.all(np.isfinite(variables_matriz)):
    print("Encontrados valores no finitos, reemplazando con 0...")
    variables_matriz[~np.isfinite(variables_matriz)] = 0

# Guardamos los nombres de los países
nombres_paises = df_limpios.columns.to_numpy()

print(f"Dimensiones de 'variables_matriz': {variables_matriz.shape}")
print(f"Longitud de 'tiempo_0': {len(tiempo_0)}")

# ==========================================================================
# 3. VISUALIZACIÓN GENERAL (Equivalente a matplot)
# ==========================================================================

# Muestra de 30 países
np.random.seed(123)
indices_muestra = np.random.choice(variables_matriz.shape[1], 30, replace=False)

# Gráfico en bruto (matplot)
plt.figure(figsize=(10, 6))
# plt.plot(x, y) grafica las columnas de y vs x, igual que matplot
plt.plot(tiempo_0, variables_matriz[:, indices_muestra])
plt.xlabel("Año")
plt.ylabel("Emisiones Totales")
plt.title("Curvas de Emisión (Muestra de 30 países aleatorios)")
plt.show()

# ==========================================================================
# 4. ANÁLISIS FDA (Equivalente a fda con scikit-fda)
# ==========================================================================

# --- Actualización variables (ya están listas) ---
tiempo = tiempo_0
variables = variables_matriz
# (nombres_paises ya está definido)

# --- Crear FDataGrid (el objeto de datos funcionales en skfda) ---
# NOTA: skfda espera los datos como (n_muestras, n_puntos_de_tiempo)
# Por lo tanto, debemos TRANSPONER nuestra matriz 'variables'.
fd_grid = skfda.FDataGrid(
    data_matrix=variables.T, # Transponer!
    grid_points=tiempo
)

# --- Creación de las bases (Equivalente a create.basis) ---
rango_tiempo = (tiempo.min(), tiempo.max())
nbasis_sugerido = 30 

# Equivalente a create.bspline.basis
base_bspline = skfda.representation.basis.BSpline(
    domain_range=rango_tiempo, 
    n_basis=nbasis_sugerido
)

# Equivalente a create.fourier.basis
# n_basis=30 es 1 (const) + 14 pares sin/cos + 1 sin extra
base_fourier = skfda.representation.basis.Fourier(
    domain_range=rango_tiempo,
    n_basis=nbasis_sugerido
)

# --- Ajustar coeficientes (Equivalente a Data2fd) ---
# .to_basis() suaviza los datos de la grilla para ajustarlos a la base
repre_bspline = fd_grid.to_basis(base_bspline)
repre_fourier = fd_grid.to_basis(base_fourier)

# --- Gráficas (Equivalente a plot.fd) ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

repre_bspline.plot(axes=axes[0])
axes[0].set_title("Gráfica DF - B-Spline")
axes[0].set_xlabel("Año")
axes[0].set_ylabel("Emisiones")

repre_fourier.plot(axes=axes[1])
axes[1].set_title("Gráfica DF - Fourier")
axes[1].set_xlabel("Año")
axes[1].set_ylabel("Emisiones")

plt.tight_layout()
plt.show()

# --- Evaluación del Error (MSE) ---
# Equivalente a eval.fd()
y_pred_bspline = repre_bspline.evaluate(tiempo).T # Transponer de (n_puntos, n_curvas) a (n_curvas, n_puntos)
y_pred_fourier = repre_fourier.evaluate(tiempo).T

mse_bspline = np.mean((variables - y_pred_bspline)**2)
mse_fourier = np.mean((variables - y_pred_fourier)**2)

print(f"MSE B-Spline: {mse_bspline:.2f}")
print(f"MSE Fourier: {mse_fourier:.2f}")

# --- Media y Varianza Funcional (B-SPLINE) ---
# Equivalente a mean.fd()
mean_fd = repre_bspline.mean()
# Equivalente a var.fd() -> skfda usa .cov()
cov_fd = repre_bspline.cov()

# --- Gráfico de curvas, media y varianza ---
grid_time = np.linspace(rango_tiempo[0], rango_tiempo[1], 200)

# fd_values_bspline <- eval.fd(grid_time, repre.bspline)
# Transponemos para que sea (n_puntos, n_curvas)
fd_values_bspline = repre_bspline.evaluate(grid_time).T 

# mean_values_bspline <- eval.fd(grid_time, mean_fd)
mean_values_bspline = mean_fd.evaluate(grid_time).flatten() # .flatten() para hacerlo un vector 1D

# var_values_bspline <- diag(eval.bifd(...))
var_values_bspline = cov_fd.evaluate(grid_time, grid_time).diagonal()
# Asegurarse de que no haya varianzas negativas diminutas por errores numéricos
var_values_bspline[var_values_bspline < 0] = 0
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
plt.ylabel("Emisiones") # Etiqueta corregida
plt.title("B-Spline: Curvas de Países, Media y Varianza")
plt.show()

# ==========================================================================
# 5. FPCA (Análisis de Componentes Principales Funcionales)
# ==========================================================================

# Equivalente a pca.fd()
nharm_max = nbasis_sugerido - 1
fpca_bspline = skfda.preprocessing.dim_reduction.FPCA(
    n_components=nharm_max
)
# Ajustamos sobre el objeto suavizado
fpca_bspline.fit(repre_bspline)

# --- Varianza acumulada ---
cumvar_bspline = np.cumsum(fpca_bspline.explained_variance_ratio_)
# k95_bspline <- which(...)[1]
try:
    k95_bspline = np.where(cumvar_bspline >= 0.95)[0][0] + 1 # +1 por índice 0
    print(f"B-Spline necesita {k95_bspline} componentes para ≥95% de varianza")
except IndexError:
    print("No se alcanzó el 95% de varianza (esto es raro)")

# --- Graficar autofunciones ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 'plot.pca.fd(..., harm=1)'
fpca_bspline.components_[0].plot(axes=axes[0], color='blue')
axes[0].set_title("B-Spline - Autofunción CP1")
axes[0].set_xlabel("Año")

# 'plot.pca.fd(..., harm=2)'
fpca_bspline.components_[1].plot(axes=axes[1], color='red')
axes[1].set_title("B-Spline - Autofunción CP2")
axes[1].set_xlabel("Año")

plt.show()

# ==========================================================================
# 6. CLUSTERING FUNCIONAL (Equivalente a mclust)
# ==========================================================================

# --- Obtener Scores (Equivalente a fpca_bspline$scores) ---
scores = fpca_bspline.transform(repre_bspline)
scores_bspline = scores[:, 0:2] # Tomamos los 2 primeros

# --- Ajustar modelo (Equivalente a Mclust(G=2:6)) ---
# Probamos de 2 a 6 clusters y elegimos el mejor (menor BIC),
# que es lo que 'Mclust' hace automáticamente.
n_components_range = range(2, 7) # R's 2:6
models = [
    GaussianMixture(n, covariance_type='full', random_state=123).fit(scores_bspline)
    for n in n_components_range
]
bics = [m.bic(scores_bspline) for m in models]

# El mejor modelo es el que tiene el BIC (Bayesian Info Criterion) más bajo
best_model_index = np.argmin(bics)
mclust_result = models[best_model_index]
n_clusters_optimo = n_components_range[best_model_index]
cluster_asignado = mclust_result.predict(scores_bspline)

print(f"El número óptimo de clusters (por BIC) es: {n_clusters_optimo}")

# --- Graficar los clusters (Equivalente a plot(mclust_result, ...)) ---
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

# --- Visualización Final: Curvas por Cluster ---
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
# Usamos el 'fd_values_bspline' (matriz de n_puntos x n_curvas) de antes
colors = plt.cm.rainbow(np.linspace(0, 1, n_clusters_optimo))

for k in range(n_clusters_optimo):
    # Encontrar los índices de las curvas en este cluster
    indices_cluster = np.where(cluster_asignado == k)[0]
    
    # Graficar todas las curvas del cluster
    plt.plot(grid_time, fd_values_bspline[:, indices_cluster], color=colors[k], alpha=0.2)

# --- Agregar la media funcional de cada clúster ---
for k in range(n_clusters_optimo):
    indices_cluster = np.where(cluster_asignado == k)[0]
    
    # mean_cluster <- mean.fd(repre.bspline[indices_cluster])
    mean_cluster = repre_bspline[indices_cluster].mean()
    mean_values = mean_cluster.evaluate(grid_time).flatten()
    
    plt.plot(grid_time, mean_values, color=colors[k], linewidth=4, linestyle='--', label=f'Media Clúster {k}')

plt.xlabel("Año")
plt.ylabel("Emisiones") # Etiqueta corregida
plt.title("Curvas de Emisión Coloreadas por Clúster")
plt.legend(loc="upper left")
plt.show()


# --- Gráficas separadas por clúster ---
print("\n--- Generando gráficos separados por clúster ---")
for k in range(n_clusters_optimo):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    indices_cluster = np.where(cluster_asignado == k)[0]
    
    # 'curvas_cluster <- eval.fd(grid_time, fd_cluster)'
    curvas_cluster = fd_values_bspline[:, indices_cluster]
    
    # 'mean_cluster <- mean.fd(fd_cluster)'
    mean_cluster = repre_bspline[indices_cluster].mean()
    mean_values = mean_cluster.evaluate(grid_time).flatten()

    # Gráfico del clúster k
    ax.plot(grid_time, curvas_cluster, color=colors[k], alpha=0.4)
    
    # Agregamos la curva media
    ax.plot(grid_time, mean_values, color=colors[k], linewidth=4, label=f'Media Clúster {k}')
    
    ax.set_xlabel("Año")
    ax.set_ylabel("Emisiones")
    ax.set_title(f"Curvas de Emisión - Clúster {k} ({len(indices_cluster)} países)")
    ax.legend(loc="upper left")
    plt.show()

print("\n--- Traducción a Python completada ---")
library(fda)         # Para representacion Funcional
library(dplyr)       # Para manipulación de datos
library(readr)       # Para leer archivos CSV (en lugar de readxl)
library(tidyr)       # Para pivotar los datos (pivot_wider)
library(stringr)     # Para limpiar nombres de países
library(ggplot2)     # Para gráficos avanzados de clustering
library(mclust)      # Para clustering avanzado


head.matrix(huella_carbono_PowerQuery)




# Limpieza datos
datos_filtrados <- huella_carbono_PowerQuery %>%
  mutate(Country = str_trim(str_replace_all(Country, "\\n", ""))) %>%
  # Nos quedamos con filas que SÍ tienen un código ISO y NO son "Global"
  filter(!is.na(`ISO 3166-1 alpha-3`), Country != "Global")

# Pivotear datos
datos_wide <- datos_filtrados %>%
  pivot_wider(
    id_cols = "Year",         # Las filas serán los años
    names_from = "Country",   # Las columnas serán los países
    values_from = "Total"     # Los valores que analizaremos
  )

# Reemplazar NA con 0
datos_wide[is.na(datos_wide)] <- 0

# VARIABLES
#tiempo <- datos_wide$Year
#tiempo

 

datos_limpios <- datos_wide %>% 
  filter(Year >= 1950, Year <= 2024) # Filtramos solo los años válidos

tiempo_0 <- datos_limpios$Year
tiempo_0
print(paste("Clase de dato 'tiempo':", class(tiempo_0)))


#str(datos_limpios[, -1])



# Nuestras 'variables' (curvas). Cada COLUMNA es ahora un país, y se converitiran los datos chr en numericos
datos_numericos <- datos_limpios %>%
  mutate(across(-Year, ~ {
    # Si es una lista, toma el primer valor numérico
    if (is.list(.x)) sapply(.x, function(v) as.numeric(v[1]))
    else as.numeric(.x)
  }))

# Verificar si hubo NAs introducidos
sum(is.na(datos_numericos))



# Crear matriz numérica final
variables_matriz <- as.matrix(datos_numericos[, -1])
variables_matriz



#variables_matriz_log <- log(variables_matriz + 1) # CAMBIAR A log(...) si los datos están muy sesgados

print(paste("Clase de variables:", class(variables_matriz)))

# Guardamos los nombres de los países para el final
nombres_paises <- colnames(variables_matriz)
nombres_paises



# Asegurarnos que es numérico
storage.mode(variables_matriz) <- "double"
print(paste("Dimensiones de variables_matriz:", paste(dim(variables_matriz), collapse = " x ")))
print(paste("Clase de variables_matriz:", class(variables_matriz)))

# Guardamos los nombres de país correctamente 
nombres_paises <- colnames(variables_matriz)
length(nombres_paises)



summary(is.finite(as.matrix(variables_matriz)))
sum(!is.finite(as.matrix(variables_matriz)))
variables_matriz[!is.finite(variables_matriz)] <- 0 #NA, NaN o Inf en la matriz. Se deben limpiar

variables_matriz <- apply(variables_matriz, 2, function(col) {
  col[!is.finite(col)] <- mean(col[is.finite(col)], na.rm = TRUE)
  col
})






# VISUALICACIÓN GENERAL

# muestra de 30 países.
#set.seed(123)
indices_muestra <- sample(ncol(variables_matriz), 30)
indices_muestra

#GRÁFICO EN BRUTO
matplot(tiempo_0, variables_matriz[, indices_muestra], type = "l", lty = 1, lwd = 1.2,
        col = rainbow(length(indices_muestra)),
        xlab = "Año", ylab = "Emisiones Totales",
        main = "Curvas de Emisión (Muestra de 30 países aleatorios)")


# TOMEMOS 10 PAISES

indices_muestra_1 <- sample(ncol(variables_matriz), 10)
indices_muestra_1
matplot(tiempo_0, variables_matriz[, indices_muestra_1], type = "l", lty = 1, lwd = 1.2,
        col = rainbow(length(indices_muestra)),
        xlab = "Año", ylab = "Emisiones Totales",
        main = "Curvas de Emisión (Muestra de 10 países aleatorios)")

# TOMEMOS 5 PAISES

indices_muestra_2 <- sample(ncol(variables_matriz), 5)
indices_muestra_2
matplot(tiempo_0, variables_matriz[, indices_muestra_2], type = "l", lty = 1, lwd = 1.2,
        col = rainbow(length(indices_muestra)),
        xlab = "Año", ylab = "Emisiones Totales",
        main = "Curvas de Emisión (Muestra de 5 países aleatorios)")

# TOMEMOS 1 PAIS
indices_muestra_3 <- sample(ncol(variables_matriz), 1)
indices_muestra_3

nombre_pais <- nombres_paises[indices_muestra_3]

matplot(tiempo_0, variables_matriz[, indices_muestra_3], type = "l", lty = 1, lwd = 1.2,
        col = rainbow(length(indices_muestra)),
        xlab = "Año", ylab = "Emisiones Totales",
        main = paste("Curva de Emisión de:", nombre_pais))


###########################################
###########################################
###########################################
###########################################
###########################################
###########################################


###########################################
###########################################
###########################################
###########################################
###########################################
###########################################



# ==========================================================================
# Análisis FDA
# ==========================================================================


# actualización variables para el FDA
tiempo <- tiempo_0

variables <- variables_matriz  # ------------------------------- ANÁLISIS SIN LOGARITMO
nombres_paises <- colnames(variables) # Guardamos los nombres para el clustering

#is.numeric(variables)


# Creación de la base B-SPLINE y FOURIER


# número de bases B-spline  --- 30
# número de bases Fourier --- 30
rango_tiempo <- c(min(tiempo), max(tiempo))
nbasis_sugerido <- 30  # PROBAR CON 30


# B-spline
base.bspline <- create.bspline.basis(rangeval = rango_tiempo, nbasis = nbasis_sugerido)

# Fourier
base.fourier <- create.fourier.basis(rangeval = rango_tiempo, nbasis = nbasis_sugerido)

# Exponencial
#tasas_crecimiento <- seq(0, 0.05, length.out = nbasis_sugerido) 
#base.exponencial <- create.exponential.basis(rangeval = rango_tiempo, ratevec = tasas)

# polygonal
#base.polygonal <- create.polygonal.basis(rangeval = rango_tiempo, argvals = tiempo)



# Ajustar coeficientes en las bases funcionales
repre.bspline <- Data2fd(argvals = tiempo, y = variables, basisobj = base.bspline)

repre.fourier <- Data2fd(argvals = tiempo, y = variables, basisobj = base.fourier)



#repre.polygonal <- Data2fd(argvals = tiempo, y = variables, basisobj = base.polygonal)






# Gráficas
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(repre.bspline, main = "Gráfica B-Spline",
     xlab = "Año", ylab = "Emisiones")
plot(repre.fourier, main = "Gráfica Fourier",
     xlab = "Año", ylab = "Emisiones")




# Es mejor B-spline que Fourier debido a que los datos no muestran periodicidad
# Idea: probar con otras bases, o cambiar el número de bases B-spline





# Evaluación del Error (MSE)
y_pred_bspline <- eval.fd(tiempo, repre.bspline)
mse_bspline <- mean((variables - y_pred_bspline)^2)

y_pred_fourier <- eval.fd(tiempo, repre.fourier)
mse_fourier <- mean((variables - y_pred_fourier)^2)



cat("MSE B-Spline:", mse_bspline, "\n")
cat("MSE Fourier:", mse_fourier, "\n")
#cat("MSE polygonal:", mse_polygonal, "\n")


# B-Spline tiene un MSE mucho menor







# Media y Varianza Funcional (B-SPLINE)

# Nos enfocaremos en B-Spline

mean_fd <- mean.fd(repre.bspline)

var_fd  <- var.fd(repre.bspline)


# Gráfico de curvas, media y varianza
grid_time <- seq(min(tiempo), max(tiempo), length.out = 200)
fd_values_bspline <- eval.fd(grid_time, repre.bspline)
mean_values_bspline <- eval.fd(grid_time, mean_fd)
cov_bspline <- eval.bifd(grid_time, grid_time, var_fd)
var_values_bspline <- diag(cov_bspline)



par(mfrow = c(1,1))
matplot(grid_time, fd_values_bspline, type = "l", lty = 1, col = "lightgray",
        xlab = "Año", ylab = "Log(Emisiones)",
        main = "B-Spline: Curvas de Países, Media y Varianza")
lines(grid_time, mean_values_bspline, col = "blue", lwd = 3)
lines(grid_time, mean_values_bspline + sqrt(var_values_bspline), col = "red", lwd = 2, lty = 2)
lines(grid_time, mean_values_bspline - sqrt(var_values_bspline), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Países", "Media", "Media ± 1 SD"),
       col = c("lightgray", "blue", "red"), lty = c(1,1,2), lwd = c(1,3,2))



# Para visualizar mejor los datos, puedo escalar 'variables' usando la función log(...)








# Análisis de Componentes Principales Funcionales (FPCA)


# Usamos B-Spline
nharm_max <- nbasis_sugerido - 1
fpca_bspline <- pca.fd(repre.bspline, nharm = nharm_max)

# Varianza acumulada
cumvar_bspline <- cumsum(fpca_bspline$varprop)
k95_bspline <- which(cumvar_bspline >= 0.95)[1]

cat("B-Spline necesita", k95_bspline, "componentes para ≥95% de varianza\n")



# Graficar autofunciones - CON DOS O 3 COMPONENTES DEBERÍA SER SUFICIENTE
#
par(mfrow = c(2,2), mar = c(4,4,2,1))
plot.pca.fd(fpca_bspline, harm = 1, nx = 200, lwd = 2, col = "blue")
#title("B-Spline - Autofunción CP1")
plot.pca.fd(fpca_bspline, harm = 2, nx = 200, lwd = 2, col = "red")
#title("B-Spline - Autofunción CP2")


















# Clustering Funcional (mclust) ---- Gaussian mixture model clustering

# Usaremos los scores de los 2 primeros componentes si no se usa el logaritmo
scores_bspline <- fpca_bspline$scores[, 1:2]
colnames(scores_bspline) <- c("CP1", "CP2")



# Ajustar modelo de clustering
mclust_result <- Mclust(scores_bspline, G = 2:6) # Buscamos de 2 a 6 clusters
summary(mclust_result) # Ver el número óptimo de clusters (G)


# Graficar los clusters en el espacio CP1-CP2
plot(mclust_result, what = "classification",
     main = "Clusters de Países (mclust)")



# Obtener el número de clusters óptimo encontrado
G_optimo <- mclust_result$G


# leyenda
legend("topleft",                                   # Posición (puedes cambiar a "topright", "bottom", etc.)
       legend = paste("Cluster", 1:G_optimo),       # Texto: "Cluster 1", "Cluster 2"...
       col = mclust.options("classPlotColors")[1:G_optimo], # Mismos colores que mclust
       pch = mclust.options("classPlotSymbols")[1:G_optimo],# Mismos símbolos que mclust
       title = "Leyenda",
       bg = "white",                                # Fondo blanco para que se lea bien
       cex = 0.8)                                   # Tamaño del texto





# Visualización Final: Curvas por Cluster
cluster_asignado <- mclust_result$classification

# Crear un dataframe para ver qué país quedó en qué cluster
df_clusters_paises <- data.frame(
  Pais = nombres_paises,
  Cluster = cluster_asignado,
  Score_CP1 = scores_bspline[, "CP1"],
  Score_CP2 = scores_bspline[, "CP2"]
)

# Ver los países de cada cluster
print("Países por cluster:")
print(split(df_clusters_paises$Pais, df_clusters_paises$Cluster))

# Graficar las curvas originales coloreadas por el cluster
par(mfrow = c(1,1))
matplot(grid_time, fd_values_bspline, type = "l", lty = 1,
        col = cluster_asignado, lwd = 1.5,
        xlab = "Año", ylab = "Log(Emisiones)",
        main = "Curvas de Emisión Coloreadas por Clúster (mclust)")

# Agregar la media funcional de cada clúster
for (k in 1:mclust_result$G) {
  mean_cluster <- mean.fd(repre.bspline[cluster_asignado == k])
  lines(grid_time, eval.fd(grid_time, mean_cluster),
        col = k, lwd = 4, lty = 2)
}
legend("topleft", legend = c(paste("Clúster", 1:mclust_result$G), "Medias"),
       col = c(1:mclust_result$G, 1), lty = c(rep(1, mclust_result$G), 2),
       lwd = c(rep(1.5, mclust_result$G), 4))









# Gráficas separadas por clúster
#mclust_result$G

# Definimos colores base
colores_clusters <- rainbow(mclust_result$G)

# Creamos una figura por clúster
for (k in 1:mclust_result$G) {
  indices_cluster <- which(cluster_asignado == k)
  
  # Subconjunto de curvas y cálculo de la media funcional del clúster
  fd_cluster <- repre.bspline[indices_cluster]
  curvas_cluster <- eval.fd(grid_time, fd_cluster)
  mean_cluster <- mean.fd(fd_cluster)
  mean_values <- eval.fd(grid_time, mean_cluster)
  
  # Gráfico del clúster k
  par(mfrow = c(1,1))
  matplot(grid_time, curvas_cluster, type = "l", lty = 1,
          col = adjustcolor(colores_clusters[k], alpha.f = 0.4),
          xlab = "Año", ylab = "Log(Emisiones)",
          main = paste("Curvas de Emisión - Clúster", k))
  
  # Agregamos la curva media
  lines(grid_time, mean_values, col = colores_clusters[k], lwd = 3)
  
  # Agregamos leyenda con los nombres de los países en el clúster
  legend("topleft",
         legend = c(paste("Clúster", k, "(", length(indices_cluster), "países)"),
                    "Media funcional"),
         col = c(adjustcolor(colores_clusters[k], alpha.f = 0.4), colores_clusters[k]),
         lty = c(1,1), lwd = c(1,3))
  
  # Mostrar nombres de países del cluster en consola
  cat("\n====================\n")
  cat("Cluster", k, "contiene los países:\n")
  cat(paste(df_clusters_paises$Pais[indices_cluster], collapse = ", "), "\n")
  cat("====================\n")
}















# ---- k-means -- #

# si no se ocupa el logaritmo
scores_bspline <- fpca_bspline$scores[, 1:2]
colnames(scores_bspline) <- c("CP1", "CP2")



# NUMERO OPTIMO CLUSTERS
set.seed(123)
wss <- numeric(10) # Vector para guardar la "Suma de Cuadrados Intra-cluster"

for (k in 1:10) {
  wss[k] <- kmeans(scores_bspline, centers = k, nstart = 25)$tot.withinss
}

# Graficamos el codo: k óptimo se encuentra cuando la línea deja de bajar bruscamente
par(mfrow = c(1,1))
plot(1:10, wss, type = "b", pch = 19, col = "blue",
     xlab = "Número de Clusters (k)", 
     ylab = "WSS (Error Intra-cluster)",
     main = "Método del Codo")

# k-optimo.... ¿se puede cambiar?
k_optimo <- 4

set.seed(123)
kmeans_result <- kmeans(scores_bspline, 
                        centers = k_optimo, 
                        iter.max = 50, 
                        nstart = 25,
                        algorithm = "Hartigan-Wong")


# Resumen estadístico real de kmeans
print(kmeans_result) 

# Guardar el cluster en un dataframe
df_scores <- data.frame(
  Pais = nombres_paises,
  CP1 = scores_bspline[, 1],
  CP2 = scores_bspline[, 2],
  Cluster = as.factor(kmeans_result$cluster) # Convertir a factor para colores
)

# GRAFICAR CLUSTERS

# Centroides (sin logaritmo)
centros <- as.data.frame(kmeans_result$centers)
colnames(centros) <- c("CP1", "CP2")


centros$Cluster <- as.factor(1:k_optimo)






ggplot(df_scores, aes(x = CP1, y = CP2, color = Cluster)) +
  # Los puntos (países)
  geom_point(alpha = 0.6, size = 2) +
  # Los centroides (X negras grandes)
  geom_point(data = centros, aes(x = CP1, y = CP2), 
             shape = 3, size = 5, stroke = 2, color = "black") +
  # Elipses (Forma del grupo)
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  labs(title = paste("K-Means Clustering (k =", k_optimo, ")"),
       subtitle = "Agrupación según patrones de emisión (CP1 y CP2)",
       x = "Componente Principal 1", 
       y = "Componente Principal 2") +
  theme_minimal()






# GRAFICAR CURVAS DE CADA CLUSTER

# ¿qué significa cada grupo?
colores <- c("red", "blue", "green3", "orange", "purple")[1:k_optimo]


dev.new(width = 12, height = 6) # Abre una ventana nueva de 12x6 pulgadas
par(mfrow = c(1, k_optimo), mar = c(2, 2, 3, 1))


     # GRAFICO POR CADA CLUSTER
for (i in 1:k_optimo) {
  indices <- which(kmeans_result$cluster == i)
  
  # Graficamos
  matplot(tiempo, variables[, indices], type = "l", lty = 1,
          col = adjustcolor(colores[i], alpha.f = 0.2), 
          xlab = "", ylab = "", # Quitamos etiquetas ejes para ahorrar espacio
          main = paste("Cluster", i, "(", length(indices), "Países)"),
          ylim = range(variables))
  
  # Agregamos la media
  media_cluster <- mean.fd(repre.bspline[indices])
  lines(media_cluster, col = "black", lwd = 3)
  
  
  # Países por cluster
  cat("\n\n======================================================\n")
  cat(paste(">>> CLUSTER", i, "- Total de países:", length(indices), "\n"))
  cat("======================================================\n")
  
  nombres_cluster <- nombres_paises[indices]
  cat(paste(nombres_cluster, collapse = ", "))
  cat("\n")
  
}



######## FIN CLUSTER K-MEANS



















# -- DBSCAN -- #

library(dbscan)



# PREPARAR DATOS --- scores del FPCA (CP1 y CP2)

# Si no se ocupa el logaritmo
scores_dbscan <- fpca_bspline$scores[, 1:2]
colnames(scores_dbscan) <- c("CP1", "CP2")




# 'EPS' ÓPTIMO (k-NN Distance Plot)

# DBSCAN necesita dos parámetros:
# 1. minPts: Cuántos puntos mínimos hacen un "barrio denso". 
#    Regla general: dimensiones + 1. Aquí tenemos 2 dim, así que minPts = 3 o 4.
#    Usaremos minPts = 4.

# 2. eps: El radio de ese barrio.


# Graficamos la distancia a los vecinos más cercanos para encontrar el "codo"
par(mfrow = c(1,1))
kNNdistplot(scores_dbscan, k = 4)

#punto donde la curva se dispara verticalmente (el "codo").
abline(h = 20, col = "red", lty = 3)


# Define tu eps basado en el gráfico anterior
mi_eps <- 10000000000   # <--- ¡CAMBIA ESTE VALOR SEGÚN TU GRÁFICO kNN!
mi_minPts <- 20



# EJECUTAR DBSCAN
set.seed(123)
dbscan_result <- dbscan(scores_dbscan, eps = mi_eps, minPts = mi_minPts)

# Ver resultados
print(dbscan_result)

# NOTA IMPORTANTE: 
# - El Cluster '0' no es un grupo real. Son los OUTLIERS (Ruido).
# - Los Clusters 1, 2, etc. son los grupos densos encontrados.

# Guardar resultados en dataframe
df_dbscan <- data.frame(
  Pais = nombres_paises,
  CP1 = scores_dbscan[, 1],
  CP2 = scores_dbscan[, 2],
  Cluster = as.factor(dbscan_result$cluster)
)




# GRAFICAR SCORES (VISUALIZACIÓN DE DENSIDAD)

ggplot(df_dbscan, aes(x = CP1, y = CP2, color = Cluster, label = Pais)) +
  geom_point(size = 2, alpha = 0.7) +
  # Usamos 'hullplot' conceptualmente con ggplot:
  # El Cluster 0 (Ruido) lo pintamos de negro o gris para destacar que son atípicos
  scale_color_manual(values = c("0" = "black", "1" = "red", "2" = "green3", "3" = "blue", "4" = "orange")) +
  labs(title = paste("DBSCAN Clustering (eps =", mi_eps, ", minPts =", mi_minPts, ")"),
       subtitle = "Cluster 0 = Ruido / Outliers (Países con comportamiento único)",
       x = "Score CP1", y = "Score CP2",
       color = "Cluster (0=Ruido)") +
  theme_minimal()




# VER QUÉ PAÍSES SON "RUIDO" (Cluster 0)
cat("\n--- PAÍSES DETECTADOS COMO OUTLIERS (RUIDO) ---\n")
outliers <- df_dbscan$Pais[df_dbscan$Cluster == 0]
print(paste(outliers, collapse = ", "))





# GRAFICAR CURVAS POR GRUPO (Separando Ruido de Grupos)

# Obtenemos los clusters únicos encontrados (incluyendo el 0)
grupos_encontrados <- sort(unique(dbscan_result$cluster))

par(mfrow = c(1, length(grupos_encontrados)), mar = c(2,2,3,1))

for (g in grupos_encontrados) {
  indices <- which(dbscan_result$cluster == g)
  
  # Título especial para el grupo 0
  titulo <- if(g == 0) paste("RUIDO / OUTLIERS (", length(indices), ")") else paste("Cluster", g, "(", length(indices), ")")
  color_grupo <- if(g == 0) "black" else rainbow(length(grupos_encontrados)-1)[g]
  
  # Graficar curvas
  matplot(tiempo, variables[, indices], type = "l", lty = 1,
          col = adjustcolor(color_grupo, alpha.f = 0.3),
          xlab = "", ylab = "",
          main = titulo,
          ylim = range(variables))
  
  # Solo dibujamos la media si NO es el grupo de ruido (el promedio de ruido no tiene sentido)
  if (g != 0) {
    media_g <- mean.fd(repre.bspline[indices])
    lines(media_g, col = "black", lwd = 3)
  }
}
par(mfrow=c(1,1)) # Reset














# -- ST-DBSCAN -- #

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)


# OBTENER COORDENADAS GEOGRÁFICAS (LAT/LON)

# Descargamos los datos del mundo
mundo <- ne_countries(scale = "medium", returnclass = "sf")


# Vamos a extraer los ISOs únicos que estamos usando en el análisis.
isos_analisis <- datos_filtrados %>% 
  select(Country, `ISO 3166-1 alpha-3`) %>% 
  distinct()

# Cruzamos con el mapa para obtener latitud y longitud del centroide de cada país
datos_geo <- mundo %>%
  select(iso_a3, geometry) %>%
  # Calculamos el centroide (el punto medio del país)
  mutate(
    lon = sf::st_coordinates(sf::st_centroid(geometry))[,1],
    lat = sf::st_coordinates(sf::st_centroid(geometry))[,2]
  ) %>%
  as.data.frame() %>%
  select(iso_a3, lon, lat)

# Unimos las coordenadas con tus datos del FPCA
df_fpca <- data.frame(
  Pais = nombres_paises,
  CP1 = scores_bspline[, 1],
  CP2 = scores_bspline[, 2]
)

# Hacemos el join usando el nombre del país o necesitamos el ISO en df_fpca
df_fpca <- df_fpca %>%
  left_join(isos_analisis, by = c("Pais" = "Country")) %>%
  left_join(datos_geo, by = c("ISO 3166-1 alpha-3" = "iso_a3"))

# Limpieza: Eliminar países que no encontraron coordenadas (si los hay)
df_st <- na.omit(df_fpca)

print(head(df_st))





# PREPARAR LA MATRIZ ESPACIO-TEMPORAL

# Aquí está el truco del ST-DBSCAN.
# Creamos una matriz con 4 dimensiones: CP1, CP2 (Tiempo/Forma) y Lat, Lon (Espacio).

matriz_st <- df_st %>% select(CP1, CP2, lon, lat)

# IMPORTANTE: ESCALADO (Scaling)
# Las coordenadas van de -180 a 180. Los Scores pueden ir de -500 a 500.
# Si no escalamos, la dimensión más grande dominará el clustering.
# scale() pone todo en unidades de "desviaciones estándar".

matriz_st_scaled <- scale(matriz_st)




# ENCONTRAR EPS ÓPTIMO (kNN Dist Plot)


par(mfrow = c(1,1))
kNNdistplot(matriz_st_scaled, k = 5)
abline(h = 1.5, col = "red", lty = 2) # Línea de referencia (ajusta según lo que veas)

# Mira el gráfico. Donde la curva sube rápido es tu epsilon.
# Supongamos que es 1.5 (al estar escalado, eps suele estar entre 0.5 y 2.0)
eps_st <- 1.5 
minPts_st <- 5






# EJECUTAR ST-DBSCAN
set.seed(123)
st_dbscan_res <- dbscan(matriz_st_scaled, eps = eps_st, minPts = minPts_st)

print(st_dbscan_res)

# Guardamos el cluster en el dataframe
df_st$ClusterST <- as.factor(st_dbscan_res$cluster)







# VISUALIZACIÓN EN UN MAPA MUNDIAL

# Cargamos el mapa base de nuevo para ggplot
world_map <- map_data("world")


ggplot() +
  # Dibujar el mapa base en gris claro
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region),
           fill = "lightgray", color = "white") +
  
  # Dibujar los puntos de nuestros países coloreados por Cluster ST
  geom_point(data = df_st, aes(x = lon, y = lat, color = ClusterST), 
             size = 3, alpha = 0.8) +
  
  # Colores (El Cluster 0 es negro/ruido)
  scale_color_manual(values = c("0" = "black", "1" = "red", "2" = "blue", 
                                "3" = "green", "4" = "orange", "5" = "purple")) +
  
  labs(title = "ST-DBSCAN: Clusters Funcionales y Geográficos",
       subtitle = "Agrupa países por cercanía geográfica Y similitud en emisiones",
       color = "Cluster ST (0=Ruido)") +
  theme_minimal()







# COMPARACIÓN

# Veamos las curvas promedio de estos grupos Espacio-Temporales

grupos <- sort(unique(st_dbscan_res$cluster))
colores <- c("black", "red", "blue", "green", "orange", "purple")

par(mfrow = c(1, length(grupos)), mar = c(2,2,3,1))

for (g in grupos) {
  # Filtramos nombres de países en este grupo
  paises_grupo <- df_st$Pais[df_st$ClusterST == g]
  
  # Encontramos los índices en la matriz original 'variables'
  indices <- which(nombres_paises %in% paises_grupo)
  
  if(length(indices) > 0) {
    titulo <- if(g==0) "Ruido / Dispersos" else paste("Cluster ST", g)
    
    matplot(tiempo, variables[, indices], type = "l", lty = 1,
            col = adjustcolor(colores[g+1], alpha.f = 0.4),
            main = titulo, ylim = range(variables), ylab="", xlab="")
    
    # Media (solo si no es ruido)
    if(g != 0) {
      media <- mean.fd(repre.bspline[indices])
      lines(media, col = "black", lwd = 3)
    }
  }
}




















#
# OBTENER COORDENADAS GEOGRÁFICAS (LAT/LON)
# Descargamos los datos del mundo
mundo <- ne_countries(scale = "medium", returnclass = "sf")

# Tu dataset original 'datos_filtrados' tiene la columna "ISO 3166-1 alpha-3".
# Vamos a extraer los ISOs únicos que estamos usando en el análisis.
isos_analisis <- datos_filtrados %>% 
  select(Country, `ISO 3166-1 alpha-3`) %>% 
  distinct()

# Cruzamos con el mapa para obtener latitud y longitud del centroide de cada país
datos_geo <- mundo %>%
  select(iso_a3, geometry) %>%
  # Calculamos el centroide (el punto medio del país)
  mutate(
    lon = sf::st_coordinates(sf::st_centroid(geometry))[,1],
    lat = sf::st_coordinates(sf::st_centroid(geometry))[,2]
  ) %>%
  as.data.frame() %>%
  select(iso_a3, lon, lat)

# Unimos las coordenadas con tus datos del FPCA
# Asumimos que 'nombres_paises' está en el mismo orden que 'scores_bspline'
df_fpca <- data.frame(
  Pais = nombres_paises,
  CP1 = scores_bspline[, 1],
  CP2 = scores_bspline[, 2]
)

# Hacemos el join usando el nombre del país o necesitamos el ISO en df_fpca
# Para hacerlo robusto, recuperemos el ISO de tus datos originales asociados a 'nombres_paises'
df_fpca <- df_fpca %>%
  left_join(isos_analisis, by = c("Pais" = "Country")) %>%
  left_join(datos_geo, by = c("ISO 3166-1 alpha-3" = "iso_a3"))

# Limpieza: Eliminar países que no encontraron coordenadas (si los hay)
df_st <- na.omit(df_fpca)

print(head(df_st))







# PREPARAR LA MATRIZ ESPACIO-TEMPORAL
# Aquí está el truco del ST-DBSCAN.
# Creamos una matriz con 4 dimensiones: CP1, CP2 (Tiempo/Forma) y Lat, Lon (Espacio).

matriz_st <- df_st %>% select(CP1, CP2, lon, lat)

# IMPORTANTE: ESCALADO (Scaling)
# Las coordenadas van de -180 a 180. Los Scores pueden ir de -500 a 500.
# Si no escalamos, la dimensión más grande dominará el clustering.
# scale() pone todo en unidades de "desviaciones estándar".

matriz_st_scaled <- scale(matriz_st)








# ENCONTRAR EPS ÓPTIMO (kNN Dist Plot)

par(mfrow = c(1,1))
kNNdistplot(matriz_st_scaled, k = 5)
abline(h = 1.5, col = "red", lty = 2) # Línea de referencia (ajusta según lo que veas)

# Mira el gráfico. Donde la curva sube rápido es tu epsilon.
# Supongamos que es 1.5 (al estar escalado, eps suele estar entre 0.5 y 2.0)
eps_st <- 1.5 
minPts_st <- 5





# EJECUTAR ST-DBSCAN
set.seed(123)
st_dbscan_res <- dbscan(matriz_st_scaled, eps = eps_st, minPts = minPts_st)

print(st_dbscan_res)

# Guardamos el cluster en el dataframe
df_st$ClusterST <- as.factor(st_dbscan_res$cluster)





# VISUALIZACIÓN EN UN MAPA MUNDIAL
# Esto es lo más valioso: ver los clusters en el mapa.

# Cargamos el mapa base de nuevo para ggplot
world_map <- map_data("world")

# Corregir algunos nombres para que coincidan con el mapa si es necesario
# Pero usaremos las coordenadas (lat/lon) que ya tenemos, así que no importa el nombre.

ggplot() +
  # Dibujar el mapa base en gris claro
  geom_map(data = world_map, map = world_map,
           aes(long, lat, map_id = region),
           fill = "lightgray", color = "white") +
  
  # Dibujar los puntos de nuestros países coloreados por Cluster ST
  geom_point(data = df_st, aes(x = lon, y = lat, color = ClusterST), 
             size = 3, alpha = 0.8) +
  
  # Colores (El Cluster 0 es negro/ruido)
  scale_color_manual(values = c("0" = "black", "1" = "red", "2" = "blue", 
                                "3" = "green", "4" = "orange", "5" = "purple")) +
  
  labs(title = "ST-DBSCAN: Clusters Funcionales y Geográficos",
       subtitle = "Agrupa países por cercanía geográfica Y similitud en emisiones",
       color = "Cluster ST (0=Ruido)") +
  theme_minimal()





#  COMPARARACIÓN

# Veamos las curvas promedio de estos grupos Espacio-Temporales

grupos <- sort(unique(st_dbscan_res$cluster))
colores <- c("black", "red", "blue", "green", "orange", "purple")

par(mfrow = c(1, length(grupos)), mar = c(2,2,3,1))

for (g in grupos) {
  # Filtramos nombres de países en este grupo
  paises_grupo <- df_st$Pais[df_st$ClusterST == g]
  
  # Encontramos los índices en la matriz original 'variables'
  indices <- which(nombres_paises %in% paises_grupo)
  
  if(length(indices) > 0) {
    titulo <- if(g==0) "Ruido / Dispersos" else paste("Cluster ST", g)
    
    matplot(tiempo, variables[, indices], type = "l", lty = 1,
            col = adjustcolor(colores[g+1], alpha.f = 0.4),
            main = titulo, ylim = range(variables), ylab="", xlab="")
    
    # Media (solo si no es ruido)
    if(g != 0) {
      media <- mean.fd(repre.bspline[indices])
      lines(media, col = "black", lwd = 3)
    }
  }
}

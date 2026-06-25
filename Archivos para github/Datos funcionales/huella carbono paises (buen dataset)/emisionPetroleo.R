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
    values_from = "Oil"     # Los valores serán las emisiones totales
  )

# Reemplazar NA con 0
datos_wide[is.na(datos_wide)] <- 0

# VARIABLES
#tiempo <- datos_wide$Year
#tiempo



datos_limpios <- datos_wide %>% 
  filter(Year >= 1900, Year <= 2024) # Filtramos solo los años válidos

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
#variables <- log(variables_matriz+1) # ------------------------------- AGREGUÉ EL LOGARITMO
nombres_paises <- colnames(variables) # Guardamos los nombres para el clustering

#is.numeric(variables)


# Creación de la base B-SPLINE y FOURIER


# número de bases B-spline  --- 30
# número de bases Fourier --- 30
rango_tiempo <- c(min(tiempo), max(tiempo))
nbasis_sugerido <- 10  # PROBAR CON 30


# B-spline
base.bspline <- create.bspline.basis(rangeval = rango_tiempo, nbasis = nbasis_sugerido)

# Fourier
base.fourier <- create.fourier.basis(rangeval = rango_tiempo, nbasis = nbasis_sugerido)

# Exponencial
tasas_crecimiento <- seq(0, 0.05, length.out = nbasis_sugerido) 
base.exponencial <- create.exponential.basis(rangeval = rango_tiempo, ratevec = tasas)

# polygonal
base.polygonal <- create.polygonal.basis(rangeval = rango_tiempo, argvals = tiempo)



# Ajustar coeficientes en las bases funcionales
repre.bspline <- Data2fd(argvals = tiempo, y = variables, basisobj = base.bspline)

repre.fourier <- Data2fd(argvals = tiempo, y = variables, basisobj = base.fourier)



repre.polygonal <- Data2fd(argvals = tiempo, y = variables, basisobj = base.polygonal)






# Gráficas
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(repre.bspline, main = "Gráfica B-Spline",
     xlab = "Año", ylab = "Emisiones")
plot(repre.fourier, main = "Gráfica Fourier",
     xlab = "Año", ylab = "Emisiones")


par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(repre.exponencial, main = "Gráfica Bases Exponenciales",
     xlab = "Año", ylab = "Emisiones")
plot(repre.polygonal, main = "Gráfica poligonal",
     xlab = "Año", ylab = "Emisiones")



# Es mejor B-spline que Fourier debido a que los datos no muestran periodicidad
# Idea: probar con otras bases, o cambiar el número de bases B-spline





# Evaluación del Error (MSE)
y_pred_bspline <- eval.fd(tiempo, repre.bspline)
mse_bspline <- mean((variables - y_pred_bspline)^2)

y_pred_fourier <- eval.fd(tiempo, repre.fourier)
mse_fourier <- mean((variables - y_pred_fourier)^2)


y_pred_polygonal <- eval.fd(tiempo, repre.polygonal)
mse_polygonal <- mean((variables - y_pred_polygonal)^2)

cat("MSE B-Spline:", mse_bspline, "\n")
cat("MSE Fourier:", mse_fourier, "\n")
cat("MSE polygonal:", mse_polygonal, "\n")
# B-Spline probablemente tiene un MSE mucho menor







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
par(mfrow = c(1,2), mar = c(4,4,2,1))
plot.pca.fd(fpca_bspline, harm = 1, nx = 200, lwd = 2, col = "blue")
#title("B-Spline - Autofunción CP1")
plot.pca.fd(fpca_bspline, harm = 2, nx = 200, lwd = 2, col = "red")
#title("B-Spline - Autofunción CP2")






# --------------------------------
# PROBANDO CON LA BASE POLYGONAL
# --------------------------------
mean_fd <- mean.fd(repre.polygonal)
var_fd  <- var.fd(repre.polygonal)

# Gráfico de curvas, media y varianza
grid_time <- seq(min(tiempo), max(tiempo), length.out = 200)
fd_values_polygonal <- eval.fd(grid_time, repre.polygonal)
mean_values_polygonal <- eval.fd(grid_time, mean_fd)
cov_polygonal <- eval.bifd(grid_time, grid_time, var_fd)
var_values_polygonal <- diag(cov_polygonal)

mean_fd

par(mfrow = c(1,1))
matplot(grid_time, fd_values_polygonal, type = "l", lty = 1, col = "lightgray",
        xlab = "Año", ylab = "Log(Emisiones)",
        main = "B-Spline: Curvas de Países, Media y Varianza")
lines(grid_time, mean_values_polygonal, col = "blue", lwd = 3)
lines(grid_time, mean_values_polygonal + sqrt(var_values_polygonal), col = "red", lwd = 2, lty = 2)
lines(grid_time, mean_values_polygonal - sqrt(var_values_polygonal), col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Países", "Media", "Media ± 1 SD"),
       col = c("lightgray", "blue", "red"), lty = c(1,1,2), lwd = c(1,3,2))



nharm_max <- nbasis_sugerido - 1
fpca_polygonal <- pca.fd(repre.polygonal, nharm = nharm_max)

# Varianza acumulada
cumvar_polygonal <- cumsum(fpca_polygonal$varprop)
k95_polygonal <- which(cumvar_polygonal >= 0.95)[1]

cat("Polygonal necesita", k95_polygonal, "componentes para ≥95% de varianza\n")

# _________________________________________________
# FIN PRUEBA CON BASE POLYGONAL --------
#__________________________________________________








# Graficar autofunciones 
#
par(mfrow = c(1,2), mar = c(4,4,2,1))
plot.pca.fd(fpca_polygonal, harm = 1, nx = 200, lwd = 2, col = "blue")
#title("B-Spline - Autofunción CP1")
plot.pca.fd(fpca_polygonal, harm = 2, nx = 200, lwd = 2, col = "red")
#title("B-Spline - Autofunción CP2")















# Clustering Funcional (mclust)
# Usaremos los scores de los 2 primeros componentes
scores_bspline <- fpca_bspline$scores[, 1:2]
colnames(scores_bspline) <- c("CP1", "CP2")


# Ajustar modelo de clustering
mclust_result <- Mclust(scores_bspline, G = 2:6) # Buscamos de 2 a 6 clusters
summary(mclust_result) # Ver el número óptimo de clusters (G)


# Graficar los clusters en el espacio CP1-CP2
plot(mclust_result, what = "classification",
     main = "Clusters de Países (mclust)")










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






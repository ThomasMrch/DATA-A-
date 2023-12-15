# Importation des librairies
library(ggplot2)
library(gridExtra)
library(car)
library(MASS)
library(PerformanceAnalytics)
library(EnvStats)
library(stats)
library(gstat)
library(data.table)
library(fields)
library(matrixStats)
library(dplyr)
library(RColorBrewer)
theme_set(theme_classic())

# Importation des données 

prediction_grid <-read.csv('prediction_grid.csv')
radio_dep <- read.csv('Radionuclide_deposition_ChNPP.csv')

head(prediction_grid)
head(radio_dep)

# Pour enlever les valeurs manquantes de Cs 137 et faire une transfo log10
radio_dep <- radio_dep[complete.cases(radio_dep$Cs.137_Bqm2), ]
radio_dep <- radio_dep %>% mutate(log_Cs.137 = log10(Cs.137_Bqm2))
radio_dep <- radio_dep %>% mutate(log_Cs.134 = log10(Cs.134_Bqm2))
radio_dep <- radio_dep %>% mutate(log_distance = log10(Distance_from_ChNPP_km))

set_cokrigeage <- radio_dep
set_cokrigeage <- set_cokrigeage[complete.cases(set_cokrigeage$Cs.134_Bqm2), ]

# Définition des intervalles
intervals <- seq(10, 360, by = 10)
intervals_rad = intervals*pi

# Création d'une liste pour stocker les données extraites et les paramètres de régression
extracted_data <- list()
regression_params <- list()
intercept_means <- numeric()

# Première passe pour calculer la moyenne des ordonnées à l'origine
for (interval in intervals) {
  # Gestion spéciale pour l'angle 10
  if (interval == 10) {
    subset_data <- radio_dep[(radio_dep$Angle_degree >= 360 - 50 & radio_dep$Angle_degree <= interval + 50) |
                               (radio_dep$Angle_degree >= interval - 50 & radio_dep$Angle_degree <= interval + 50) |
                               (radio_dep$Angle_degree >= interval - 50 + 360 & radio_dep$Angle_degree <= interval + 50 + 360), ]
  } else {
    # Gestion spéciale pour l'angle 360
    if (interval == 360) {
      subset_data <- radio_dep[(radio_dep$Angle_degree >= interval - 50 & radio_dep$Angle_degree <= interval + 50) |
                                 (radio_dep$Angle_degree >= interval - 50 - 360 & radio_dep$Angle_degree <= interval + 50 - 360) |
                                 (radio_dep$Angle_degree >= interval - 50 + 360 & radio_dep$Angle_degree <= interval + 50), ]
    } else {
      # Sélection des données pour l'intervalle spécifié
      subset_data <- radio_dep[(radio_dep$Angle_degree >= interval - 50 & radio_dep$Angle_degree <= interval + 50) |
                                 (radio_dep$Angle_degree >= interval - 50 + 350 & radio_dep$Angle_degree <= interval + 50 - 360), ]
    }
  }
  
  # Stockage des données extraites dans la liste
  extracted_data[[as.character(interval)]] <- subset_data
  
  # Régression linéaire avec Distance_from_ChNPP_km
  lm_model <- lm(log_Cs.137 ~ log_distance, data = subset_data)
  
  # Stockage des paramètres de régression dans la liste
  regression_params[[as.character(interval)]] <- coef(lm_model)
  
  # Calcul de la moyenne des ordonnées à l'origine
  intercepts <- sapply(regression_params, function(params) params[1])
  intercept_mean <- mean(intercepts)
  intercept_means <- c(intercept_means, intercept_mean)
}

# Deuxième passe pour ajuster les modèles avec ordonnée à l'origine fixe
regression_params_fixed_intercept <- list()
for (interval in intervals) {
  # Extraire les données
  subset_data <- extracted_data[[as.character(interval)]]
  # Ajuster le modèle linéaire avec ordonnée à l'origine fixe
  lm_model_fixed_intercept <- lm(log_Cs.137 ~ 0+ log_distance, offset = rep(intercept_mean, nrow(subset_data)), data = subset_data)
  #lm_model_fixed_intercept <- lm(I(log_Cs.137-intercept_mean) ~ 0 + Distance_from_ChNPP_km, data = subset_data)
  
  
  # Stocker les paramètres de régression dans la liste
  regression_params_fixed_intercept[[as.character(interval)]] <- coef(lm_model_fixed_intercept)
  #print(coef(lm_model_fixed_intercept))
}


# Extraction des pentes des régressions
slopes_sans_offset <- sapply(regression_params, function(params) params[2])
slopes_avec_offset <-sapply(regression_params_fixed_intercept, function(params) params[1])


derniere_pente_avec_offset <- slopes_avec_offset[length(slopes_avec_offset)]
slopes_avec_offset <- c(derniere_pente_avec_offset, slopes_avec_offset)

intervals_bis <- seq(0, 360, by = 10)
intervals_bis_rad = intervals_bis*pi/180

# Ajout d'un point supplémentaire au début des données pour contraindre la spline à commencer en 0
angles_interpolation <- seq(0, 360, by = 0.5) * pi /180
slopes_interpolation_avec_offset <- spline(intervals_bis_rad, slopes_avec_offset, n = length(angles_interpolation))$y


#------- Si on ajuste correctement la spline avec 3 x 36 valeurs

#faire 3 x 36 points donc on va avoir un cycle qui se rép
#graphique de visualisation de 3 periodes

# Répliquer les données deux fois de plus pour obtenir trois périodes
intervals_triplicated <- c(intervals_bis, intervals_bis + 360, intervals_bis + 720)
slopes_triplicated <- c(slopes_avec_offset, slopes_avec_offset, slopes_avec_offset)

# Adapter les angles en radians pour les trois périodes
intervals_triplicated_rad <- intervals_triplicated * pi / 180

# Interpolation des données pour obtenir une fonction continue sur les trois périodes
angles_interpolation_triplicated <- seq(0, 1080, by = 0.5) * pi / 180  # Trois périodes complètes
slopes_interpolation_triplicated <- spline(intervals_triplicated_rad, slopes_triplicated, n = length(angles_interpolation_triplicated))$y

#isolement de la periode du milieu 

# Sélection des données pour la période de 2π à 4π radians
angles_2pi_to_4pi <- angles_interpolation_triplicated[seq(720, 1440)]  # Sélection des angles de 2π à 4π
slopes_2pi_to_4pi <- slopes_interpolation_triplicated[seq(720, 1440)]  # Sélection des pentes correspondantes

# Ajuster les angles pour correspondre à un vrai cercle (de 0 à 2π)
angles_2pi_to_4pi_adjusted <- angles_2pi_to_4pi - 2 * pi  # Ajustement des angles de 2π à 4π à une plage de 0 à 2π


adjusted_angle_rad <- prediction_grid$Angle_rad + 2 * pi

# Utilisez la fonction spline avec les points spécifiques (Xout)
spline_result_modified <- spline(x = intervals_triplicated_rad, y = slopes_triplicated, xout = adjusted_angle_rad)

# Ajoutez les valeurs interpolées modifiées à prediction_grid
prediction_grid$Interpolated_Slope_Modified <- spline_result_modified$y
prediction_grid$Esperance_Modified <- log10(prediction_grid$Distance_from_ChNPP_km) * prediction_grid$Interpolated_Slope_Modified + intercept_mean

#---------------Variogramme--------------

log_Cs137.gstat <- gstat(formula = log_Cs.137~1, data = radio_dep, locations = ~longitude+latitude)

distance_matrix <- dist(radio_dep[, c("longitude", "latitude")], method = "euclidean") #dist max pour statuer la demi moitié
hmax <- max(distance_matrix)
var_log_Cs137 = var(radio_dep$log_Cs.137,na.rm=T)

print(paste("Maximum distance between two points : ",round(hmax,3))) #pour arrondir
print(paste("Variance of head :",round(var_log_Cs137,3)))

log_Cs137.vario <- variogram(log_Cs137.gstat)
log_Cs137.vario #faire en sorte que chaque np soit supérieur à 30, pour chaque distance

ggplot(log_Cs137.vario) + 
  geom_point(aes(x = dist,y=gamma))+
  labs(x='distance',y='semivariance') +
  ggtitle("Variogram log Cs 137") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(aes(yintercept = var_log_Cs137, linetype="Variance")) +
  scale_x_continuous(limits = c(0, max(log_Cs137.vario$dist))) +
  scale_y_continuous(limits = c(0, var_log_Cs137+0.35)) ##Pour forcer l'échelle à 0 en y mais jusqu'à la valeur de la variance +1000 pour que la variance s'affiche

# Faire la même chose mais avec les résidus
#en gros on recalcule l'espérance selon nos points de log_Cs137 dans radio_dep

adjusted_angle_radio_dep <- radio_dep$Angle_rad + 2 * pi
spline_radio_dep_log <- spline(x = intervals_triplicated_rad, y = slopes_triplicated, xout = adjusted_angle_radio_dep)
radio_dep$spline_radio_dep_log <- spline_radio_dep_log$y
radio_dep$Esperance_log <- radio_dep$log_distance * radio_dep$spline_radio_dep_log + intercept_mean

#Calcul des résidus = E[x] = obs - résidus ==> résidus = obs - E[x]
radio_dep$Log_residus <- radio_dep$log_Cs.137 - radio_dep$Esperance_log


log_Cs137_var.res.gstat <- gstat(formula = Log_residus~1, data = radio_dep, locations = ~longitude+latitude)

distance_matrix <- dist(radio_dep[, c("longitude", "latitude")], method = "euclidean") #dist max pour statuer la demi moitié
hmax <- max(distance_matrix)
var_res_log_Cs137 = var(radio_dep$Log_residus,na.rm=T)

print(paste("Maximum distance between two points : ",round(hmax,3))) #pour arrondir
print(paste("Variance of residuals head :",round(var_res_log_Cs137,3)))

res_log_Cs137.vario <- variogram(log_Cs137_var.res.gstat)
res_log_Cs137.vario #faire en sorte que chaque np soit supérieur à 30, pour chaque distance

ggplot(res_log_Cs137.vario) + 
  geom_point(aes(x = dist,y=gamma))+
  labs(x='distance',y='semivariance') +
  ggtitle("Variogram log Cs-137 résidus") + #Becquerels par metre carré
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(aes(yintercept = var_res_log_Cs137, linetype="Variance")) +
  scale_x_continuous(limits = c(0, max(res_log_Cs137.vario$dist))) +
  scale_y_continuous(limits = c(0, var_res_log_Cs137+0.1)) ##Pour forcer l'échelle à 0 en y mais jusqu'à la valeur de la variance +1000 pour que la variance s'affiche


#####______ fitted variogram ________

vario.res_log_Cs137 <- variogram(log_Cs137_var.res.gstat)
model.res_log_c137 <- vgm(psill=0.05, model='Sph', range = 0.3, nugget = 0.12) # parameter values estimated from experimental variogram

plot(vario.res_log_Cs137,model.res_log_c137)

model.res_log_Cs137.fit <- fit.variogram(vario.res_log_Cs137, model = model.res_log_c137)
print(model.res_log_Cs137.fit)

plot(vario.res_log_Cs137,model.res_log_Cs137.fit)


#####___________Krigeage_____________

# Kriging for Cs 137
Log_Cs137.krig <- krige(formula = log_Cs.137 ~ 1, data = radio_dep, locations = ~longitude + latitude, newdata = prediction_grid, model = model.res_log_Cs137.fit)

# Create a plot for Cesium 137 kriging prediction
ggplot() + 
  geom_tile(data = Log_Cs137.krig, aes(x = longitude, y = latitude, fill = var1.pred)) +
  geom_point(data = radio_dep, aes(x = longitude, y = latitude, color = ""), shape = 18, size = 2) +
  scale_color_manual("", values = "black") +
  scale_fill_gradientn(name = "", colors = c('royalblue', 'turquoise', 'yellow', 'red')) +
  theme(legend.key = element_rect(fill = "transparent", color = NA)) +
  labs(x = "Longitude (°W)", y = "Latitude (°N)", col = NULL) + 
  ggtitle("Prédictions des concentrations en Log 137 - Cs [log(Bq/m²)]") + #log Becquerels par metre carré
  theme(plot.title = element_text(hjust = 0.5))

ggplot() + 
  geom_tile(data = Log_Cs137.krig, aes(x = longitude, y = latitude, fill = var1.var)) +
  geom_point(data = radio_dep, aes(x = longitude, y = latitude, color = ""), shape = 18, size = 2) +
  scale_color_manual("", values = "black") +
  scale_fill_gradientn(name = "", 
                       colors = c('royalblue', 'turquoise', 'yellow', 'red')) +
  theme(legend.key = element_rect(fill = "transparent", color = NA)) +
  labs(x = "Longitude (°W)", y = "Latitude (°N)", col = NULL) + 
  ggtitle("Variance des prédictions de Log 137 - Cs") + 
  theme(plot.title = element_text(hjust = 0.5))


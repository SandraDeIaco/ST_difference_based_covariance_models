##### Load packages #####

library(dplyr)
library(lattice)
library(RColorBrewer)
library(readxl)
library(nloptr)

##### Load the empirical spatio-temporal covariance ##### 

load("cova_st.RData")

##### 3D plot of empirical spatio-temporal covariance ##### 

wireframe(cova_st[,3] ~ cova_st[,1]*cova_st[,2],
          data = cova_st,
          scales = list(arrows = FALSE, cex = 0.55),
          drape = TRUE, colorkey = TRUE,
          xlab = list("spatial lags", rot = -56),
          ylab = list("temporal lags", rot = 19),
          zlab = list(label = "", rot = 90),
          screen = list(z = -65, x = -65),
          perspective = TRUE,
          zlim = c(-0.2, 1.01)
          )


##### 2D plot of empirical spatial covariance ##### 

cova_hs <- cova_st[cova_st$timelag == 0, ]
plot(x = cova_hs$spacelag, y = cova_hs$cova)


##### 2D plot of empirical temporal covariance ##### 
cova_ht <- cova_st[cova_st$spacelag == 0, ]
plot(x = cova_ht$timelag, y = cova_ht$cova)


##### Spatio-temporal covariance model with optimized parameters ##### 

###### Step 1: optimization of the spatial marginal parameters ######
source("Step1_C_space_bis.R", echo=TRUE)

#Save the estimated parameters
A1 <- res_sp$solution[1]
B1 <- res_sp$solution[2]
c1 <- res_sp$solution[3]
c2 <- res_sp$solution[4]
c11 <- res_sp$solution[5]
c21 <- res_sp$solution[6]


# #test on the spatial parameters
# c1/c2 < (A1/B1)^0.5
# c11/c21 < (A1/B1)

###### Step 2: optimization of the temporal marginal parameters ######
source("Step2_C_time_bis.R", echo=TRUE)

#Save the estimated parameters
A1 <- res_tm$solution[1]
B1 <- res_tm$solution[2]
a1 <- res_tm$solution[3]
a2 <- res_tm$solution[4]

 
# #test on the temporal parameters
# a1/a2
# (A1/B1)^2
# a1/a2 < (A1/B1)^2
# 
# a2/a1
# c21/c11
# (a2/a1) < (c21/c11)


###### Step 3: optimization of the spatio-temporal parameters ######
source("Step3_C_spacetime_bis.R", echo=TRUE)

#Save the estimated parameters
A1 <- res_st$solution[1]
B1 <- res_st$solution[2]
c1 <- res_st$solution[3]
c2 <- res_st$solution[4]
c11 <- res_st$solution[5]
c21 <- res_st$solution[6]
a1 <- res_st$solution[7]
a2 <- res_st$solution[8]

# #test on the spatio-temporal parameters
# 1 <= (c1/c2) 
# (c1/c2) < (A1/B1)^(1/2)
# 
# 
# 1 < (c11/c21)^2
# (c11/c21)^2 < (a1/a2) * (c11/c21)
# (a1/a2) * (c11/c21) < (A1/B1)^2


###### Computation of the spatio-temporal theoretical covariance, according to the estimated parameters ###### 
cova_st_th <- st_model(h = cova_st$spacelag, u = cova_st$timelag, A1 = A1, B1 = B1, c1 = c1, c2 = c2, 
                       c11 = c11, c21 = c21, a1 = a1, a2 = a2, 
                       alpha1 = 0.1, alpha2 = 0.9)


cova_st_th <- data.frame(spacelag = cova_st$spacelag, timelag = cova_st$timelag,
                         covateo = cova_st_th)

###### 3D plot of spatio-temporal theoretical covariance
wireframe(cova_st_th[,3] ~ cova_st_th[,1] * cova_st_th[,2],
          data = cova_st_th,
          scales = list(arrows = FALSE, cex = 0.55),
          drape = TRUE, colorkey = TRUE,
          xlab = list("spatial lags", rot = -56),
          ylab = list("temporal lags", rot = 19),
          zlab = list(label = "", rot = 90),
          screen = list(z = -65, x = -65),
          perspective = TRUE,
          zlim = c(-0.25, 1.01)
          )

###### Computation of error metrics ######
cova_st_th_emp <- merge(cova_st_th, cova_st, by = c("spacelag", "timelag"))
cova_st_th_emp$abs_sc <- abs(cova_st_th_emp$covateo - cova_st_th_emp$cova)
cova_st_th_emp$sq_sc <- (cova_st_th_emp$covateo - cova_st_th_emp$cova)^2

MAE <- mean(cova_st_th_emp$abs_sc) #mean absolute error
RMSE <- sqrt(mean(cova_st_th_emp$sq_sc)) #root mean square error


######2D plot of empirical and theoretical spatial and temporal covariances ######
h <- seq(from = 0.0, to = 180, by = 1)
u <- seq(from = 0.0, to = 9, by = 0.1)
lag_st <- expand.grid(h, u)
names(lag_st) <- c("h", "u")
cova_st_th <- st_model(h = lag_st$h, u = lag_st$u, A1 = A1, B1 = B1, c1 = c1, c2 = c2, 
                       c11 = c11, c21 = c21, a1 = a1, a2 = a2, 
                       alpha1 = 0.1, alpha2 = 0.9)
cova_st_th <- data.frame(spacelag = lag_st$h, timelag = lag_st$u,
                         covateo = cova_st_th)

sp_th <- subset(cova_st_th, timelag == 0)
plot(x = cova_hs$spacelag, y = cova_hs$cova, ylim = c(-0.2, 1))
lines(h, sp_th$covateo, col = "red", lwd = 2)


tm_th <- subset(cova_st_th, spacelag == 0)
plot(x = cova_ht$timelag, y = cova_ht$cova)
lines(u, tm_th$covateo, col = "red", lwd = 2)

# =========================================================
# STEP 3 - Spatio-temporal parameters fitting
#
# A1, B1, a1 and a2 are inherited from Step2_C_time_is.R
#
# c1, c2, c11 and c21 are inherited from Step1_C_space.R
# estimated parameters: A1, B1, c1, c2, c11, c21, a1, a2.
# =========================================================


# Spatial domain dimension
#d <- 2 


# Data
h_data    <- cova_st$spacelag
u_data    <- cova_st$timelag
cova_data <- cova_st$cova


# =======================================================
#Theoretical covariance model - spatio-temporal - eq (63) 
# =======================================================

st_model <- function(h, u, A1 = res_tm$solution[1], B1 = res_tm$solution[2], 
                     c1 = res_sp$solution[3], c2 = res_sp$solution[4],
                     c11 = res_sp$solution[5], c21 = res_sp$solution[6], 
                     a1 = res_tm$solution[3], a2 = res_tm$solution[4], 
                     alpha1 = 0.1, alpha2 = 0.9) {
  
  alpha1 * (A1 * exp(-c1 * h) - B1 * exp(-c2 * h)) + 
    alpha2 * ((A1 /(1 + a1 * abs(u)^2)) * exp(-((c11 * h)^2)/(1 + a1 * abs(u)^2)) - (B1 /(1 + a2 * abs(u)^2)) * exp(-((c21 * h)^2)/(1 + a2 * abs(u)^2)))
}


# =========================================================
# Objective function
# Sum of Squared Errors
# =========================================================
loss_function_st <- function(params, h, u, cova_emp) {
  A1 <- params[1]
  B1 <- params[2]
  c1 <-  params[3]
  c2 <-  params[4]
  c11 <- params[5]
  c21 <- params[6]
  a1 <- params[7]
  a2 <- params[8]
  
  # Theoretical spatio-temporal covariance model - eq (61) 
  cova_th_st <- st_model(h, u, A1, B1, c1, c2, c11, c21, a1, a2, alpha1 = 0.1, alpha2 = 0.9)
  
  # Error Sum of Squares
  sse <- sum((cova_emp - cova_th_st)^2)
  
  sse 
}



# ==============================================================================
# Inequality constraints
#
# nloptr convention:
# g_i(x) <= 0
#
# Constraints:
# 1 <= (c1/c2) < (A1/B1)^(1/2)
# 1 < (c11/c21)^d < (a1/a2) * (c11/c21)^(d-1) < (A1/B1)^2
# ==============================================================================
constraints_st <- function(params) {
  eps <- 1e-4
  A1 <- params[1]
  B1 <- params[2] + eps
  c1 <-  params[3]
  c2 <- params[4]
  c11 <- params[5]
  c21 <- params[6]
  a1 <- params[7]
  a2 <- params[8]
  
  c(
  
    # -----------------------------------------
    # Constraint 1:
    # 1 <= c1/c2
    # (c1/c2) < (A1/B1)^(1/2)
    # -----------------------------------------
    1 - (c1/c2),
    (c1/c2) - (A1/B1)^(1/2),
    
    
    # -----------------------------------------
    # Constraint 2:
    # 1 < (c11/c21)^d
    # (c11/c21)^d < (a1/a2) * (c11/c21)^(d-1)
    # (a1/a2) * (c11/c21)^(d-1) < (A1/B1)^2
    # -----------------------------------------
    1 - (c11/c21)^2,
    (c11/c21)^2 - (a1/a2) * (c11/c21),
    (a1/a2) * (c11/c21) - (A1/B1)^2
  
  )

}


# =========================================================
# Starting values
# c1, c2, c11, c21 inherited from Step1_C_space.R
# A1, B1, a1, a2 inherited from Step2_C_time.R
# =========================================================

start_vals_st <- c(A1 = res_tm$solution[1], B1 = res_tm$solution[2],
                   c1 = res_sp$solution[3], c2 = res_sp$solution[4],
                   c11 = res_sp$solution[5], c21 = res_sp$solution[6],
                   a1 = res_tm$solution[3], a2 = res_tm$solution[4])


#Main function for parameters optimization (nloptr::nloptr)
res_st <- nloptr::nloptr(
  x0 = start_vals_st,
  eval_f = function(params) loss_function_st(params, 
                                             h = h_data, 
                                             u = u_data,
                                             cova_emp = cova_data),
  
  eval_g_ineq = constraints_st,
  lb = rep(1e-6, 8),
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-8, maxeval = 10000)
)

res_st$solution

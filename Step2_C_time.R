# =========================================================
# STEP 2 - Temporal marginal parameters fitting
#
# A1 and B1 are inherited from Step1_C_space.R
# Estimated parameters: A1, B1, a1 and a2.
# =========================================================


# Data
u_data    <- cova_ht$timelag
cova_data <- cova_ht$cova


# =========================================================
# Theoretical covariance model - temporal marginal (eq. 65)
# =========================================================

tm_marg_model <- function(u, A1 = res_sp$solution[1], 
                          B1 = res_sp$solution[2], 
                          a1, a2, alpha1 = 0.1, alpha2 = 0.9) {
  alpha1 + alpha2 * ((A1 / (1 + a1 * abs(u)^2)) - (B1 / (1 + a2 * abs(u)^2 )))
  
}


# =========================================================
# Objective function
# Sum of Squared Errors
# =========================================================

loss_function_tm <- function(params, u, cova_emp) {
  A1 <- params[1]
  B1 <- params[2]
  a1 <- params[1]
  a2 <- params[2]
  
  # Theoretical covariance model - marginal in time - eq (65)
  cova_th_tm <- tm_marg_model(u, A1, B1, a1, a2, alpha1 = 0.1, alpha2 = 0.9)
  
  # Error Sum of Squares
  sse <- sum((cova_emp - cova_th_tm)^2)
  sse
}


# =========================================================
# Inequality constraints
#
# nloptr convention:
# g_i(x) <= 0
#
# Constraints:
#
# 1 <= a1/a2 < (A1/B1)^(2)
# a2/a1 < c21/c11
# =========================================================

constraints_tm <- function(params) {
  eps <- 1e-4  
  A1 <- params[1] 
  B1 <- params[2] + eps
  a1 <- params[3] - eps
  a2 <- params[4] 
  
  c(
    
    # -----------------------------------------
    # Constraint 1:
    # 1 <= a1/a2
    # -----------------------------------------
    1 - (a1 / a2),
    

    # a1/a2 <= (A1/B1)^2
    (a1 / a2) - (A1 / B1)^2,
    
    # -----------------------------------------
    # Constraint 2:
    # a2/a1 < c21/c11
    (a2 / a1) - (c21 / c11)
    
  )
}


# =========================================================
# Starting values
# A1 and B1 inherited from Step1_C_space.R
# =========================================================

start_vals_tm <- c(
  A1 = res_sp$solution[1],
  B1 = res_sp$solution[2],
  a1 = 1.5, a2 = 0.1)


# =========================================================
# Optimization
# =========================================================

res_tm <- nloptr::nloptr(
  x0 = start_vals_tm,
  eval_f = function(params) {loss_function_tm(params = params, 
                                              u = u_data, 
                                              cova_emp = cova_data)
    },
  
  eval_g_ineq = constraints_tm,
  lb = rep(1e-6, 4),
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-8, maxeval = 10000)
)


# =========================================================
# Estimated parameters
# =========================================================

res_tm$solution
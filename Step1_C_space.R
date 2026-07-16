# =========================================================
# STEP 1 - Spatial marginal parameters fitting
# =========================================================


# Spatial domain dimension
#d <- 2


# Data
h_data    <- cova_hs$spacelag
cova_data <- cova_hs$cova


# =========================================================
# Theoretical covariance model - marginal in space (eq. 64)
# =========================================================

sp_marg_model <- function(h, A1, B1, c1, c2, c11, c21, 
                          alpha1 = 0.1, alpha2 = 0.9) {
  alpha1 * (A1 * exp(-c1 * h) - B1 * exp(-c2 * h)) + 
    alpha2 * (A1 * exp(-(c11 * h)^2) - B1 * exp(-(c21 * h)^2))
}


#Loss function
loss_function_sp <- function(params, h, cova_emp) {
  A1 <- params[1]
  B1 <- params[2]
  c1 <- params[3]
  c2 <- params[4]
  c11 <- params[5]
  c21 <- params[6]
  
  # Theoretical covariance model - marginal in space - eq (64) 
  cova_th_sp <- sp_marg_model(h, A1, B1, c1, c2, c11, c21, alpha1 = 0.1, alpha2 = 0.9)
  
  # Error Sum of Squares
  sse <- sum((cova_emp - cova_th_sp)^2)
  
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
# 1 <= c1/c2 < (A1/B1)^(1/2)
#
# 1 <= c11/c21 <= (A1/B1)
# =========================================================

constraints_sp <- function(params) {
  eps <- 1e-4
  A1  <- params[1]
  B1  <- params[2] + eps
  c1  <- params[3]
  c2  <- params[4]
  c11 <- params[5]
  c21 <- params[6]
  
  c(
    
    # -----------------------------------------
    # Constraint 1
    # 1 <= c1/c2
    # -----------------------------------------
    1 - (c1 / c2),
    

    # c1/c2 < (A1/B1)^(1/2)
    (c1 / c2) - (A1 / B1)^(1/2),
    
    # -----------------------------------------
    # Constraint 2
    # 1 <= c11/c21
    # -----------------------------------------
    1 - (c11 / c21),
    

    # c11/c21 <= A1/B1
    (c11 / c21) - (A1 / B1)
    
  )
}


# =========================================================
# Starting values
# =========================================================
start_vals_sp <- c(A1 = 1.0, B1 = 0.1, c1 = 0.1, c2 = 0.01, 
                   c11 = 0.21, c21 = 0.01)


# =========================================================
# Optimization
# =========================================================
#Main function for parameters optimization (nloptr::nloptr)
res_sp <- nloptr::nloptr(
  x0 = start_vals_sp,
  eval_f = function(params) {loss_function_sp(params = params,
                                              h = h_data, 
                                              cova_emp = cova_data)
    },
  eval_g_ineq = constraints_sp,
  lb = rep(1e-6, 6),
  opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-8, maxeval = 10000)
)

res_sp$solution
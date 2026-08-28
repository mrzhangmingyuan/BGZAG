library(stats)
library(numDeriv)
library(MASS)
library(Matrix)

set.seed(123)

# ===================
# 基础函数定义
# ===================

rSGe <- function(n, p) {
  rgeom(n, prob = p) + 1
}

generate_Z <- function(n, a) {
  sample(c(1, 0, -1), size = n, replace = TRUE, 
         prob = c(a^2, 2*a*(1-a), (1-a)^2))
}

# 模拟函数
simulate_model <- function(n = 1500, omega =3, alpha = 0.1, beta = 0.3,
                           gamma = -0.7, a = 0.6, delta = 5) {
  X <- numeric(n)
  Y <- numeric(n)
  Z <- numeric(n)
  p <- numeric(n)
  lambda <- numeric(n)
  
  # 初始化
  lambda[1] <- ifelse(1 - beta != 0, omega / (1 - beta), 1e-8)
  phi <- 2 * a^2 - 2 * a + 1
  
  for (t in 2:n) {
    # 更新lambda
    lambda[t] <- (omega + alpha * (abs(X[t - 1]) - gamma * X[t - 1])^delta + beta * lambda[t - 1]^delta)^(1/delta)
    lambda[t] <- max(lambda[t], 1e-8)
    
    # 计算p
    p[t] <- (sqrt(phi^2 + 4 * lambda[t]) - phi) / (2 * lambda[t])
    p[t] <- min(max(p[t], 1e-8), 1 - 1e-8)
    
    # 生成Y和Z
    Y[t] <- rSGe(1, p[t])
    Z[t] <- generate_Z(1, a)
    
    # 计算X
    X[t] <- Z[t] * Y[t]
  }
  
  return(list(X = X, Y = Y, Z = Z, p = p, lambda = lambda))
}

# 似然函数
log_likelihood <- function(params, X) {
  n <- length(X)
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  gamma <- params[4]
  a <- params[5]
  delta <- params[6]
  
  # 参数约束
  if (omega <= 0 || alpha < 0 || beta < 0 || 
      gamma < -1 || gamma > 1 || a <= 0 || a >= 1 || 
      delta <= 0 || delta > 10) {
    return(-Inf)
  }
  if (alpha + beta >= 1) {
    return(-Inf)
  }
  
  lambda <- numeric(n)
  p <- numeric(n)
  log_lik <- 0
  
  lambda[1] <- omega / (1 - beta)
  phi <- 2 * a^2 - 2 * a + 1
  
  for (t in 2:n) {
    lambda[t] <- (omega + alpha * (abs(X[t - 1]) - gamma * X[t - 1])^delta + 
                    beta * lambda[t - 1]^delta)^(1/delta)
    lambda[t] <- max(lambda[t], 1e-8)
    
    p[t] <- (sqrt(phi^2 + 4 * lambda[t]) - phi) / (2 * lambda[t])
    p[t] <- min(max(p[t], 1e-8), 1 - 1e-8)
    
    if (X[t] == 0) {
      log_lik <- log_lik + log(2 * a * (1 - a))
    } else {
      y_abs <- abs(X[t])
      z_sign <- sign(X[t])
      
      geom_prob <- (1 - p[t])^(y_abs - 1) * p[t]
      
      if (z_sign == 1) {
        z_prob <- a^2
      } else if (z_sign == -1) {
        z_prob <- (1 - a)^2
      } else {
        z_prob <- 2 * a * (1 - a)
      }
      
      log_lik <- log_lik + log(geom_prob) + log(z_prob)
    }
    
    if (!is.finite(log_lik)) {
      return(-Inf)
    }
  }
  return(log_lik)
}

# 参数估计函数
estimate_parameters_simple <- function(X, start_params = NULL) {
  if (is.null(start_params)) {
    start_params <- c(omega = 3, alpha = 0.1, beta = 0.3, 
                      gamma = -0.7, a = 0.6, delta = 5)
  }
  
  lower_bounds <- c(0.001, 0, 0, -0.99, 0.01, 0.01)
  upper_bounds <- c(5, 0.99, 0.99, 0.99, 0.99, 10)
  
  objective <- function(params) {
    ll <- log_likelihood(params, X)
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  res <- optim(par = start_params, fn = objective, 
               method = "L-BFGS-B", 
               lower = lower_bounds, upper = upper_bounds,
               control = list(maxit = 500))
  
  if (res$convergence != 0) return(NULL)
  
  names(res$par) <- c("omega", "alpha", "beta", "gamma", "a", "delta")
  return(list(parameters = res$par, log_likelihood = -res$value, convergence = res$convergence))
}

# ===================
# Monte Carlo 仿真
# ===================

true_params <- c(omega = 3, alpha = 0.1, beta = 0.3, gamma = -0.7, a = 0.6, delta = 5)

num_iterations <- 500   # 可改大
estimated_params_list <- list()

for (i in 1:num_iterations) {
  sim_data <- simulate_model(n = 1500, 
                             omega = true_params[1], 
                             alpha = true_params[2], 
                             beta = true_params[3], 
                             gamma = true_params[4], 
                             a = true_params[5], 
                             delta = true_params[6])
  
  res <- tryCatch({
    estimate_parameters_simple(sim_data$X)
  }, error = function(e) NULL)
  
  if (!is.null(res) && res$convergence == 0) {
    estimated_params_list[[length(estimated_params_list) + 1]] <- res$parameters
  } else {
    cat(sprintf("第 %d 次优化失败，跳过。\n", i))
  }
}

# ===================
# 结果整理与汇总
# ===================

if (length(estimated_params_list) > 0) {
  results_matrix <- do.call(rbind, estimated_params_list)
  
  means <- colMeans(results_matrix)
  sds <- apply(results_matrix, 2, sd)
  ses <- sds / sqrt(nrow(results_matrix))
  bias <- means - true_params
  rmse <- sqrt(colMeans((results_matrix - 
                           matrix(true_params, nrow(results_matrix), 
                                  ncol = length(true_params), byrow = TRUE))^2))
  
  summary_table <- data.frame(
    Parameter = names(true_params),
    True_Value = true_params,
    Mean_Estimate = round(means, 4),
    Std_Dev = round(sds, 4),
    Std_Error = round(ses, 4),
    Bias = round(bias, 4),
    RMSE = round(rmse, 4)
  )
  
  print(summary_table, row.names = FALSE)
} else {
  cat("警告: 所有优化都失败，没有结果可汇总。\n")
}

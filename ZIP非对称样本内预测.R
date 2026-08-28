library(stats)
library(numDeriv)
library(MASS)
library(readxl)

# ---- 1. 读取数据 ----
data <- read_excel("D:/万科A股.xlsx") 

X_data <- data$y  # 假设收益列名为y

# ---- 2. 对称 ZIP 分布 ----
sym_zip_pmf <- function(k, omega, lambda) {
  if (k == 0) {
    return(omega + (1 - omega) * exp(-lambda))
  } else {
    return(0.5 * (1 - omega) * dpois(abs(k), lambda))
  }
}

sym_zip_cdf <- function(k, omega, lambda) {
  if (k < 0) {
    return(0.5 * (1 - omega) * ppois(abs(k), lambda))
  } else if (k == 0) {
    return(omega + (1 - omega) * exp(-lambda) +
             0.5 * (1 - omega) * (1 - exp(-lambda)))
  } else {
    return(omega + (1 - omega) * exp(-lambda) +
             0.5 * (1 - omega) * (ppois(k, lambda) + ppois(k, lambda) - 1))
  }
}

# ---- 3. 负对数似然 ----
negative_log_likelihood <- function(params, X_data) {
  omega_zip <- params[1]  
  alpha <- params[2]
  beta <- params[3]
  gamma <- params[4]
  delta <- params[5]
  omega_garch <- params[6]  
  
  if (omega_zip <= 0 || omega_zip >= 1 || alpha < 0 || beta < 0 || beta >= 1 || 
      delta <= 0 || omega_garch <= 0) {
    return(1e10)  
  }
  
  n <- length(X_data)
  lambda <- numeric(n)
  lambda[1] <- ifelse(1 - beta != 0, omega_garch / (1 - beta), 1e-8)
  log_likelihood <- 0
  
  for (t in 2:n) {
    term1 <- abs(X_data[t - 1]) - gamma * X_data[t - 1]
    if (term1 < 0) term1 <- 1e-8  
    
    lambda[t] <- (omega_garch + alpha * term1^delta + beta * lambda[t - 1]^delta)^(1/delta)
    lambda[t] <- max(lambda[t], 1e-8)  
    
    k <- X_data[t]
    prob <- sym_zip_pmf(k, omega_zip, lambda[t])
    
    if (is.na(prob) || prob <= 0) {
      log_likelihood <- log_likelihood - 1e6  
    } else {
      log_likelihood <- log_likelihood + log(prob)
    }
  }
  
  if (!is.finite(log_likelihood)) {
    return(1e10)
  }
  
  return(-log_likelihood)  
}

# ---- 4. 参数设置 ----
init_params <- c(omega_zip = 0.3, alpha = 0.1, beta = 0.2, gamma = 0.6, delta = 0.5, omega_garch = 0.1)

# ---- 5. 优化 ----
result <- optim(
  par = init_params,
  fn = negative_log_likelihood,
  X_data = X_data,
  method = "Nelder-Mead",
  control = list(maxit = 2000),
  hessian = TRUE
)

# ---- 6. 标准误 ----
if (result$convergence == 0) {
  hessian <- result$hessian
  if (all(eigen(hessian)$values > 0)) {
    cov_matrix <- solve(hessian)
    se <- sqrt(diag(cov_matrix))
  } else {
    warning("Hessian矩阵非正定，使用广义逆近似。")
    cov_matrix <- ginv(hessian)
    se <- sqrt(abs(diag(cov_matrix)))
  }
} else {
  stop("优化未收敛！")
}# ---- 7. 样本内预测 ----
set.seed(13)
n <- length(X_data)
RMSE_values <- numeric(50)
MSE_values  <- numeric(50)
MAE_values  <- numeric(50)

# 提取估计参数
params <- result$par
omega_zip   <- params[1]
alpha       <- params[2]
beta        <- params[3]
gamma       <- params[4]
delta       <- params[5]
omega_garch <- params[6]

# 循环预测 test_size = 1 到 50
for (test_size in 1:50) {
  
  # 测试集
  test_data <- X_data[(n - test_size + 1):n]
  
  # ---- 样本内拟合 λ ----
  lambda <- numeric(n)
  lambda[1] <- ifelse(1 - beta != 0, omega_garch / (1 - beta), 1e-8)
  
  for (t in 2:n) {
    term1 <- abs(X_data[t - 1]) - gamma * X_data[t - 1]
    if (term1 < 0) term1 <- 1e-8
    lambda[t] <- (omega_garch + alpha * term1^delta + beta * lambda[t - 1]^delta)^(1/delta)
    lambda[t] <- max(lambda[t], 1e-8)
  }
  
  # ---- 预测未来 test_size 步 ----
  predictions <- numeric(test_size)
  lambda_pred <- numeric(test_size)
  
  lambda_prev <- lambda[n - test_size + 1]
  X_prev <- X_data[n - test_size + 1]
  
  for (i in 1:test_size) {
    term1 <- abs(X_prev) - gamma * X_prev
    if (term1 < 0) term1 <- 1e-8
    
    lambda_pred[i] <- (omega_garch + alpha * term1^delta + beta * lambda_prev^delta)^(1/delta)
    lambda_pred[i] <- max(lambda_pred[i], 1e-8)
    
    # ---- 生成预测值：对称 ZIP 分布 ----
    u <- runif(1)
    if (u < omega_zip) {
      k_sim <- 0
    } else {
      y_sim <- rpois(1, lambda_pred[i])
      if (y_sim == 0) {
        k_sim <- 0
      } else {
        sign_choice <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
        k_sim <- sign_choice * y_sim
      }
    }
    predictions[i] <- k_sim
    
    # 更新状态
    lambda_prev <- lambda_pred[i]
    X_prev <- predictions[i]
  }
  
  # ---- 误差指标 ----
  RMSE_values[test_size] <- sqrt(mean((test_data - predictions)^2))
  MSE_values[test_size]  <- mean((test_data - predictions)^2)
  MAE_values[test_size]  <- mean(abs(test_data - predictions))
}

# ---- 8. 输出结果 ----
selected_test_sizes <- c(5, 10, 15)
cat("RMSE values:\n", RMSE_values[selected_test_sizes], "\n")
cat("MSE values:\n",  MSE_values[selected_test_sizes], "\n")
cat("MAE values:\n",  MAE_values[selected_test_sizes], "\n")


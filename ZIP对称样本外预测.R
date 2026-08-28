library(stats)
library(numDeriv)
library(MASS)
library(readxl)

# ---- 1. 读取数据 ----
data <- read_xlsx("D:/中国银行_平衡1.xlsx") 
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
}# ---- 7. 样本外预测与误差评估 ----
set.seed(1)
n <- length(X_data)

RMSE_values <- numeric(50)
MSE_values <- numeric(50)
MAE_values <- numeric(50)

for (test_size in 1:50) {
  train_data <- X_data[1:(n - test_size)]
  test_data <- X_data[(n - test_size + 1):n]
  train_size <- length(train_data)
  
  # 初始化 lambda
  lambda_train <- numeric(train_size)
  lambda_train[1] <- ifelse(1 - result$par[3] > 1e-6, result$par[6] / (1 - result$par[3]), 1e-8)
  
  # 样本内拟合 λ
  for (t in 2:train_size) {
    term1 <- abs(train_data[t - 1]) - result$par[4] * train_data[t - 1]
    if (term1 < 0) term1 <- 1e-8
    lambda_train[t] <- (result$par[6] + result$par[2] * term1^result$par[5] +
                          result$par[3] * lambda_train[t - 1]^result$par[5])^(1/result$par[5])
    lambda_train[t] <- max(lambda_train[t], 1e-8)
  }
  
  # 样本外预测
  predictions <- numeric(test_size)
  lambda_pred <- numeric(test_size)
  
  X_prev <- train_data[train_size]
  lambda_prev <- lambda_train[train_size]
  
  for (i in 1:test_size) {
    term1 <- abs(X_prev) - result$par[4] * X_prev
    if (term1 < 0) term1 <- 1e-8
    
    lambda_pred[i] <- (result$par[6] + result$par[2] * term1^result$par[5] +
                         result$par[3] * lambda_prev^result$par[5])^(1/result$par[5])
    lambda_pred[i] <- max(lambda_pred[i], 1e-8)
    
    # 从对称 ZIP 分布中模拟预测值
    if (runif(1) < result$par[1]) {
      # 零膨胀部分
      predictions[i] <- 0
    } else {
      # Poisson 对称部分
      k <- rpois(1, lambda_pred[i])
      if (k == 0) {
        predictions[i] <- 0
      } else {
        predictions[i] <- sample(c(-k, k), 1)
      }
    }
    
    # 更新
    X_prev <- predictions[i]
    lambda_prev <- lambda_pred[i]
  }
  
  # 计算误差指标
  RMSE_values[test_size] <- sqrt(mean((test_data - predictions)^2))
  MSE_values[test_size] <- mean((test_data - predictions)^2)
  MAE_values[test_size] <- mean(abs(test_data - predictions))
}

# ---- 8. 输出结果 ----
selected_test_sizes <- c(5, 10, 15)

cat("RMSE values:\n", RMSE_values[selected_test_sizes], "\n")
cat("MSE values:\n", MSE_values[selected_test_sizes], "\n")
cat("MAE values:\n", MAE_values[selected_test_sizes], "\n")

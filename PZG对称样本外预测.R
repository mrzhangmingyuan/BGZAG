library(readxl)
library(stats)
library(numDeriv)
library(MASS)
library(ggplot2)

# 读取数据
data <- read_xlsx("D:/中国银行_平衡1.xlsx") 
X_data <- data$y  # 假设收益列名为y，且为整数值

# 定义Poisson分布的概率质量函数 (PMF)
possion_pmf <- function(k, lambda) {
  k <- abs(k) 
  if(k == 0) {
    return(exp(-lambda))
  } else {
    return(0.5 * exp(-lambda) * (lambda^k) / factorial(k))
  }
}

# 定义负对数似然函数
negative_log_likelihood <- function(params, X_data) {
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  gamma <- params[4]
  
  n <- length(X_data)
  eta <- numeric(n)
  lambda <- numeric(n)
  
  # 初始化 eta 和 lambda
  eta[1] <- ifelse(1 - beta != 0, omega / (1 - beta), 1e-8)
  lambda[1] <- (sqrt(1 + 4 * eta[1]) - 1) / 2
  
  log_likelihood <- 0
  
  for (t in 2:n) {
    # 计算 eta_t
    eta[t] <- omega + alpha * (abs(X_data[t - 1]) - gamma * X_data[t - 1])^2 + beta * eta[t - 1]
    
    # 计算 lambda_t
    lambda[t] <- (sqrt(1 + 4 * eta[t]) - 1) / 2
    
    # 计算Poisson分布的PMF
    k <- X_data[t]
    prob <- possion_pmf(k, lambda[t])
    
    if (is.na(prob) || prob <= 0) {
      log_likelihood <- log_likelihood - 1e6
    } else {
      log_likelihood <- log_likelihood + log(prob)
    }
  }
  
  return(-log_likelihood)
}

# 参数初始化和边界
init_params <- c(omega = 2, alpha = 0.2, beta = 0.3, gamma =0.8)
lower_bounds <- c(1e-8, 0, 0, -1)
upper_bounds <- c(10, 1, 1, 1)

# 执行优化
result <- optim(
  par = init_params,
  fn = negative_log_likelihood,
  X_data = X_data,
  method = "Nelder-Mead",
  control = list(maxit = 2000),
  hessian = TRUE
)

# 计算标准误
if (result$convergence == 0) {
  hessian <- result$hessian
  if (all(eigen(hessian)$values > 0)) {
    cov_matrix <- solve(hessian)
    se <- sqrt(diag(cov_matrix))
  } else {
    warning("Hessian矩阵非正定，使用广义逆近似计算标准误。")
    cov_matrix <- MASS::ginv(hessian)
    se <- sqrt(diag(cov_matrix))
  }
} else {
  stop("优化未收敛，标准误不可靠！")
}
set.seed(6)
RMSE_values <- numeric(50)
MSE_values <- numeric(50)
MAE_values <- numeric(50)
n <- length(X_data)  # 数据总长度

for (test_size in 1:50) {
  
  train_data <- X_data[1:(n - test_size)]  # 训练集
  test_data <- X_data[(n - test_size + 1):n]  # 测试集
  
  train_size <- length(train_data)  # 训练集大小
  
  # 计算训练集上的 eta 和 lambda
  eta_train <- numeric(train_size)
  lambda_train <- numeric(train_size)
  
  eta_train[1] <- ifelse(1 - result$par[3] > 1e-6, result$par[1] / (1 - result$par[3]), 1e-8)
  lambda_train[1] <- (sqrt(1 + 4 * eta_train[1]) - 1) / 2
  
  for (t in 2:train_size) {
    eta_train[t] <- result$par[1] + result$par[2] * (abs(train_data[t - 1]) - result$par[4] * train_data[t - 1])^2 + result$par[3] * eta_train[t - 1]
    lambda_train[t] <- (sqrt(1 + 4 * eta_train[t]) - 1) / 2
  }
  
  # 预测未来 test_size 期
  predictions <- numeric(test_size)
  eta_pred <- numeric(test_size)
  lambda_pred <- numeric(test_size)
  
  # 预测的初始值从训练集的最后一个 lambda 开始
  eta_prev <- eta_train[train_size]
  X_prev <- train_data[train_size]
  
  for (i in 1:test_size) {
    eta_pred[i] <- result$par[1] + result$par[2] * (abs(X_prev) - result$par[4] * X_prev)^2 + result$par[3] * eta_prev
    lambda_pred[i] <- (sqrt(1 + 4 * eta_pred[i]) - 1) / 2
    predictions[i] <- rpois(1, lambda_pred[i])
    
    # 更新用于下一步预测的 X_prev 和 lambda_prev
    X_prev <- predictions[i]
  }
  
  # 计算 RMSE
  RMSE_values[test_size] <- sqrt(mean((test_data - predictions)^2))
  # 计算 MSE
  MSE_values[test_size] <- mean((test_data - predictions)^2)
  MAE_values[test_size] <- mean(abs(test_data - predictions))
}

# 选择 test_size = 5, 10, 15 的值
selected_test_sizes <- c(5, 10, 15)

cat(" RMSE values:\n", RMSE_values[selected_test_sizes], "\n")
cat(" MSE values:\n", MSE_values[selected_test_sizes], "\n")
cat(" MAE values:\n", MAE_values[selected_test_sizes], "\n")

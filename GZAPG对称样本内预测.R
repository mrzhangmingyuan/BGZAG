library(stats)
library(numDeriv)
library(MASS)
library(readxl)  # 你需要确保已加载此库来读取Excel文件

# ---- 1. 读取真实数据 ----
data <- read_xlsx("D:/中国银行_平衡1.xlsx") 
X_data <- data$y  # 假设收益列名为y

# ---- 2. 定义几何分布相关函数 ----
# 几何分布PMF：P(Y = k) = (1-p)^(k-1) * p, k = 1, 2, 3, ...
geometric_pmf <- function(k, p) {
  if (k <= 0) return(0)
  return((1 - p)^(k - 1) * p)
}

geometric_cdf <- function(k, p) {
  if (k <= 0) return(0)
  return(1 - (1 - p)^k)
}

# ---- 3. 辅助函数 ----
generate_Z <- function(n, a) {
  sample(c(1, 0, -1), size = n, replace = TRUE, 
         prob = c(a^2, 2*a*(1-a), (1-a)^2))
}

# ---- 4. 负对数似然函数（直接使用观测值计算） ----
negative_log_likelihood <- function(params, X_data) {
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  gamma <- params[4]
  a <- params[5]
  delta <- params[6]
  
  # 参数约束检查
  if (omega <= 0 || alpha < 0 || beta < 0 || beta >= 1 || 
      a <= 0 || a >= 1 || delta <= 0) {
    return(1e10)  # 返回一个很大的值
  }
  
  n <- length(X_data)
  lambda <- numeric(n)
  p <- numeric(n)
  lambda[1] <- ifelse(1 - beta != 0, omega / (1 - beta), 1e-8)
  phi <- 2 * a^2 - 2 * a + 1
  log_likelihood <- 0
  
  for (t in 2:n) {
    # 更新lambda
    term1 <- abs(X_data[t - 1]) - gamma * X_data[t - 1]
    if (term1 < 0) term1 <- 1e-8  # 避免负值的幂运算
    
    lambda[t] <- (omega + alpha * term1^delta + beta * lambda[t - 1]^delta)^(1/delta)
    lambda[t] <- max(lambda[t], 1e-8)  # 避免lambda为负
    
    # 计算p值
    p[t] <- (sqrt(phi^2 + 4 * lambda[t]) - phi) / (2 * lambda[t])
    p[t] <- min(max(p[t], 1e-8), 1 - 1e-8)  # 确保p在(0,1)范围内
    
    # 直接计算观测值X[t]的概率
    k <- X_data[t]
    
    if (k == 0) {
      prob <- 2 * a * (1 - a)  # P(Z = 0)
    } else {
      y_abs <- abs(k)
      z_sign <- sign(k)
      
      # 几何分布的概率
      geom_prob <- geometric_pmf(y_abs, p[t])
      
      # Z的概率
      if (z_sign == 1) {
        z_prob <- a^2
      } else if (z_sign == -1) {
        z_prob <- (1 - a)^2
      } else {
        z_prob <- 2 * a * (1 - a)
      }
      
      prob <- geom_prob * z_prob
    }
    
    # 处理数值下溢问题
    if (is.na(prob) || prob <= 0) {
      log_likelihood <- log_likelihood - 1e6  # 惩罚项
    } else {
      log_likelihood <- log_likelihood + log(prob)
    }
  }
  
  # 检查是否有无穷大或 NaN
  if (!is.finite(log_likelihood)) {
    return(1e10)
  }
  
  return(-log_likelihood)  # 返回负对数似然
}

# ---- 5. 参数初始化和边界 ----
init_params <- c(omega = 0.1, alpha = 0.1, beta = 0.2, gamma = 0.7, a = 0.7, delta = 0.7)
lower_bounds <- c(0, 0, 0, -1, 0, 0)
upper_bounds <- c(10, 1, 1, 1, 1, 5)

# ---- 6. 执行优化 ----
result <- optim(
  par = init_params,
  fn = negative_log_likelihood,
  X_data = X_data,
  method = "Nelder-Mead",
  control = list(maxit = 2000),
  hessian = TRUE
)

# ---- 7. 计算标准误 ----
if (result$convergence == 0) {
  hessian <- result$hessian
  if (all(eigen(hessian)$values > 0)) {
    cov_matrix <- solve(hessian)
    se <- sqrt(diag(cov_matrix))
  } else {
    warning("Hessian矩阵非正定，使用广义逆近似计算标准误。")
    cov_matrix <- ginv(hessian)
    se <- sqrt(abs(diag(cov_matrix)))  # 使用绝对值避免负方差
  }
} else {
  stop("优化未收敛，标准误不可靠！")
}# ---- 8. 样本内预测 ----
set.seed(30)
n <- length(X_data)
RMSE_values <- numeric(50)
MSE_values <- numeric(50)
MAE_values <- numeric(50)

# 循环 test_size = 1 到 50
for (test_size in 1:50) {
  
  # 提取最后 test_size 个数据作为测试集
  test_data <- X_data[(n - test_size + 1):n]
  
  # 用估计参数重新计算 lambda 和 p
  params <- result$par
  omega <- params[1]; alpha <- params[2]; beta <- params[3]
  gamma <- params[4]; a <- params[5]; delta <- params[6]
  
  lambda <- numeric(n)
  p <- numeric(n)
  phi <- 2 * a^2 - 2 * a + 1
  
  # 初始化 lambda
  lambda[1] <- ifelse(1 - beta != 0, omega / (1 - beta), 1e-8)
  
  for (t in 2:n) {
    term1 <- abs(X_data[t - 1]) - gamma * X_data[t - 1]
    if (term1 < 0) term1 <- 1e-8
    lambda[t] <- (omega + alpha * term1^delta + beta * lambda[t - 1]^delta)^(1/delta)
    lambda[t] <- max(lambda[t], 1e-8)
    
    p[t] <- (sqrt(phi^2 + 4 * lambda[t]) - phi) / (2 * lambda[t])
    p[t] <- min(max(p[t], 1e-8), 1 - 1e-8)
  }
  
  # ---- 预测 test_size 步 ----
  predictions <- numeric(test_size)
  lambda_pred <- numeric(test_size)
  p_pred <- numeric(test_size)
  
  lambda_prev <- lambda[n - test_size + 1]
  X_prev <- X_data[n - test_size + 1]
  
  for (i in 1:test_size) {
    term1 <- abs(X_prev) - gamma * X_prev
    if (term1 < 0) term1 <- 1e-8
    lambda_pred[i] <- (omega + alpha * term1^delta + beta * lambda_prev^delta)^(1/delta)
    lambda_pred[i] <- max(lambda_pred[i], 1e-8)
    
    p_pred[i] <- (sqrt(phi^2 + 4 * lambda_pred[i]) - phi) / (2 * lambda_pred[i])
    p_pred[i] <- min(max(p_pred[i], 1e-8), 1 - 1e-8)
    
    # 随机生成预测值：先生成几何分布，再乘以符号变量 Z
    Y_sim <- rgeom(1, prob = p_pred[i]) + 1
    Z_sim <- sample(c(1, 0, -1), size = 1, prob = c(a^2, 2*a*(1-a), (1-a)^2))
    predictions[i] <- Y_sim * Z_sim
    
    # 更新用于下一步预测的状态
    lambda_prev <- lambda_pred[i]
    X_prev <- predictions[i]
  }
  
  # ---- 计算误差指标 ----
  RMSE_values[test_size] <- sqrt(mean((test_data - predictions)^2))
  MSE_values[test_size]  <- mean((test_data - predictions)^2)
  MAE_values[test_size]  <- mean(abs(test_data - predictions))
}

# ---- 9. 输出结果 ----
selected_test_sizes <- c(5, 10, 15)

cat("RMSE values:\n", RMSE_values[selected_test_sizes], "\n")
cat("MSE values:\n",  MSE_values[selected_test_sizes], "\n")
cat("MAE values:\n",  MAE_values[selected_test_sizes], "\n")


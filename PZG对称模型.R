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

# 输出参数估计结果
cat("\nPZG-GARCH(1,1) 参数估计结果（基于真实数据）：\n")
cat(sprintf("omega = %.4f (SE = %.4f)\n", result$par[1], se[1]))
cat(sprintf("alpha = %.4f (SE = %.4f)\n", result$par[2], se[2]))
cat(sprintf("beta  = %.4f (SE = %.4f)\n", result$par[3], se[3]))
cat(sprintf("gamma = %.4f (SE = %.4f)\n", result$par[4], se[4]))

# 计算 AIC 和 BIC
n <- length(X_data)  # 样本数
k <- length(result$par)  # 参数个数
log_likelihood <- -result$value  # 取负值恢复对数似然

AIC <- 2 * k - 2 * log_likelihood
BIC <- k * log(n) - 2 * log_likelihood

# 输出 AIC 和 BIC
cat(sprintf("AIC  = %.4f\n", AIC))
cat(sprintf("BIC  = %.4f\n", BIC))
cat(sprintf("对数似然函数值 (log-likelihood) = %.2f\n", log_likelihood))
RMSE_values <- numeric(50)
MSE_values <- numeric(50)

n <- length(X_data)  # 数据总长度
# 提取最后 test_size 个数据作为测试集
test_data <- X_data[1:n]

lambda <- numeric(n)
eta <- numeric(n)
eta[1] <- ifelse(1 - result$par[3] != 0, result$par[1] / (1 - result$par[3]), 1e-8)
lambda[1] <- (sqrt(1 + 4 * eta[1]) - 1) / 2

for (t in 2:n) {
  eta[t] <- result$par[1] + result$par[2] * (abs(X_data[t - 1]) - result$par[4] * X_data[t - 1])^2 + result$par[3] * eta[t - 1]
  lambda[t] <- (sqrt(1 + 4 * eta[t]) - 1) / 2
}
set.seed(1234)
predictions <- numeric(n)
eta_pred <- numeric(n)
lambda_pred <- numeric(n)

eta_prev <- eta[1]
X_prev <- X_data[1]

for (i in 1:n) {
  eta_pred[i] <- result$par[1] + result$par[2] * (abs(X_prev) - result$par[4] * X_prev)^2 + result$par[3] * eta_prev
  lambda_pred[i] <- (sqrt(1 + 4 * eta_pred[i]) - 1) / 2
  
  # 生成Poisson值
  Y <- rpois(1, lambda_pred[i])
  
  # 生成 ±1 符号
  sign_Y <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
  
  # 对称预测值
  predictions[i] <- Y * sign_Y
  
  # 更新 X_prev 和 eta_prev
  X_prev <- predictions[i]
  eta_prev <- eta_pred[i]
}



# 计算 MADE

# 计算 RMSE
RMSE <- sqrt(mean((X_data - predictions)^2))
# 计算 MSE
MSE <- mean((X_data - predictions)^2)
MAE <- mean(abs(X_data - predictions))

# 输出 RMS 和 MAE
cat(sprintf("RMS (Root Mean Square) = %.4f\n", RMSE))
cat(sprintf("MAE (Mean Absolute Error) = %.4f\n", MAE))
cat(sprintf("MSE (Mean Absolute Error) = %.4f\n", MSE))
# ---- 对称 PZG-GARCH(1,1) 预测 + RQR ----
set.seed(1234)

predictions <- numeric(n)
lambda_pred <- numeric(n)
eta_pred <- numeric(n)

eta_prev <- eta[1]
X_prev <- X_data[1]

for (i in 1:n) {
  eta_pred[i] <- result$par[1] + result$par[2] * (abs(X_prev) - result$par[4] * X_prev)^2 + result$par[3] * eta_prev
  lambda_pred[i] <- (sqrt(1 + 4 * eta_pred[i]) - 1) / 2
  
  # 生成 Poisson 绝对值
  Y <- rpois(1, lambda_pred[i])
  
  # 随机符号 ±1
  sign_Y <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
  predictions[i] <- Y * sign_Y
  
  # 更新
  X_prev <- predictions[i]
  eta_prev <- eta_pred[i]
}

# ---- 计算 PIT ----
pit_values <- numeric(n)
for (t in 1:n) {
  k <- X_data[t]
  lambda_t <- lambda[t]
  
  if (k == 0) {
    pit_values[t] <- runif(1, 0, 0.5 + 0.5 * exp(-lambda_t))  # 0 对应累积概率
  } else {
    y_abs <- abs(k)
    lower <- 0.5 + 0.5 * ppois(y_abs - 1, lambda_t)  # 对称累积概率下界
    upper <- 0.5 + 0.5 * ppois(y_abs, lambda_t)      # 对称累积概率上界
    if (k < 0) {
      # 负数部分对称处理
      lower <- 0.5 - 0.5 * ppois(y_abs, lambda_t)
      upper <- 0.5 - 0.5 * ppois(y_abs - 1, lambda_t)
    }
    pit_values[t] <- runif(1, lower, upper)
  }
}

# 防止 PIT 为 0 或 1
epsilon <- 1e-8
RQR <- qnorm(pmin(pmax(pit_values, epsilon), 1 - epsilon))

# ---- 绘图 ----
par(mfrow = c(1, 2))

# PIT 直方图
hist(pit_values, breaks = 9, freq = FALSE, col = "skyblue", border = "black",
     xlab = "PIT", ylab = "Relative Frequency", main = "(b1) PIT Histogram of PZG(1,1) Model")
abline(h = 1, col = "red", lwd = 2, lty = 2)

# RQR QQ 图
qqnorm(RQR, main = "(b2) QQ Plot of RQR (PZG(1,1))")
qqline(RQR, col = "red", lwd = 2)

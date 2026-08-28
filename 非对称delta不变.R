library(stats)
library(numDeriv)
library(MASS)
library(readxl)  # 你需要确保已加载此库来读取Excel文件

# ---- 1. 读取真实数据 ----
data <- read_excel("D:/万科A股.xlsx") 
X_data <- data$y  # 假设收益列名为y
delta = 1
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

  
  # 参数约束检查
  if (omega <= 0 || alpha < 0 || beta < 0 || beta >= 1 || 
      a <= 0 || a >= 1) {
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
init_params <- c(omega = 0.5, alpha = 0.1, beta = 0.1, gamma = -0.5, a = 0.6)
lower_bounds <- c(0, 0, 0, -1, 0)
upper_bounds <- c(10, 1, 1, 1, 1)

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
}

# ---- 8. 输出结果 ----
cat("\n几何分布模型参数估计结果（基于真实数据）：\n")
cat(sprintf("omega = %.4f (SE = %.4f)\n", result$par[1], se[1]))
cat(sprintf("alpha = %.4f (SE = %.4f)\n", result$par[2], se[2]))
cat(sprintf("beta  = %.4f (SE = %.4f)\n", result$par[3], se[3]))
cat(sprintf("gamma = %.4f (SE = %.4f)\n", result$par[4], se[4]))
cat(sprintf("a     = %.4f (SE = %.4f)\n", result$par[5], se[5]))

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

set.seed(1234)

# ---- 9. 计算 lambda 和进行预测 ----
lambda <- numeric(n)
p_values <- numeric(n)
lambda[1] <- ifelse(1 - result$par[3] != 0, result$par[1] / (1 - result$par[3]), 1e-8)
phi <- 2 * result$par[5]^2 - 2 * result$par[5] + 1

for (t in 2:n) {
  term1 <- abs(X_data[t - 1]) - result$par[4] * X_data[t - 1]
  if (term1 < 0) term1 <- 1e-8
  
  lambda[t] <- (result$par[1] + result$par[2] * term1^delta + 
                  result$par[3] * lambda[t - 1]^delta)^(1/delta)
  lambda[t] <- max(lambda[t], 1e-8)  # 避免负数
  
  p_values[t] <- (sqrt(phi^2 + 4 * lambda[t]) - phi) / (2 * lambda[t])
  p_values[t] <- min(max(p_values[t], 1e-8), 1 - 1e-8)
}

# 预测未来值
predictions <- numeric(n)
lambda_pred <- numeric(n)
lambda_prev <- lambda[1]  # 预测起点的 lambda
X_prev <- X_data[1]  # 预测起点的 X

for (i in 1:n) {
  term1 <- abs(X_prev) - result$par[4] * X_prev
  if (term1 < 0) term1 <- 1e-8
  
  lambda_pred[i] <- (result$par[1] + result$par[2] * term1^delta + 
                       result$par[3] * lambda_prev^delta)^(1/delta)
  lambda_pred[i] <- max(lambda_pred[i], 1e-8)  # 避免负数
  
  p_pred <- (sqrt(phi^2 + 4 * lambda_pred[i]) - phi) / (2 * lambda_pred[i])
  p_pred <- min(max(p_pred, 1e-8), 1 - 1e-8)
  
  # 生成一个基于几何分布的预测值
  Y_pred <- rgeom(1, prob = p_pred) + 1  # 几何分布预测
  Z_pred <- sample(c(1, 0, -1), size = 1, prob = c(result$par[5]^2, 
                                                   2*result$par[5]*(1-result$par[5]), (1-result$par[5])^2))
  predictions[i] <- Y_pred * Z_pred
  
  # 更新用于下一步预测的 X_prev 和 lambda_prev
  X_prev <- predictions[i]
  lambda_prev <- lambda_pred[i]
}

# ---- 10. 计算 RMS 和 MAE ----
RMS <- sqrt(mean((X_data - predictions)^2))
MAE <- mean(abs(X_data - predictions))
MSE <- mean((X_data - predictions)^2)

# 输出 RMS 和 MAE
cat(sprintf("RMS (Root Mean Square) = %.4f\n", RMS))
cat(sprintf("MAE (Mean Absolute Error) = %.4f\n", MAE))
cat(sprintf("MSE (Mean Square Error) = %.4f\n", MSE))

# ---- 计算 lambda_hat 和 p_hat（用于 PIT 和 RQR） ----
lambda_hat <- numeric(n)
p_hat <- numeric(n)
lambda_hat[1] <- ifelse(1 - result$par[3] != 0, result$par[1] / (1 - result$par[3]), 1e-8)
a_hat <- result$par[5]

for (t in 2:n) {
  term1 <- abs(X_data[t - 1]) - result$par[4] * X_data[t - 1]
  if (term1 < 0) term1 <- 1e-8
  
  lambda_hat[t] <- (result$par[1] + result$par[2] * term1^delta + 
                      result$par[3] * lambda_hat[t - 1]^delta)^(1/delta)
  lambda_hat[t] <- max(lambda_hat[t], 1e-8)
  
  p_hat[t] <- (sqrt((2*a_hat^2 - 2*a_hat + 1)^2 + 4 * lambda_hat[t]) - (2*a_hat^2 - 2*a_hat + 1)) / (2 * lambda_hat[t])
  p_hat[t] <- min(max(p_hat[t], 1e-8), 1 - 1e-8)
}

# ---- 计算 PIT 值 ----
pit_values <- numeric(n)
for (t in 1:n) {
  k <- X_data[t]
  if (k == 0) {
    pit_values[t] <- runif(1, 0, 1)
  } else {
    y_abs <- abs(k)
    lower <- geometric_cdf(y_abs - 1, p_hat[t])
    upper <- geometric_cdf(y_abs, p_hat[t])
    pit_values[t] <- runif(1, lower, upper)
  }
}
pit_values_nonzero <- pit_values[X_data != 0]

epsilon <- 1e-8
RQR <- qnorm(pmin(pmax(pit_values, epsilon), 1 - epsilon))

# ---- 绘图 ----
par(mfrow = c(1, 2))

# PIT直方图
hist(pit_values_nonzero, breaks = 9, freq = FALSE, col = "skyblue", border = "black",
     xlab = "PIT", ylab = "Relative Frequency", 
     main = "(A1) PIT Histogram of GZAPG(1,1) Model")
abline(h = 1, col = "red", lwd = 2, lty = 2)

# RQR QQ图
qqnorm(RQR, main = "(A2) QQ Plot of RQR (GZAPG(1,1))")
qqline(RQR, col = "red", lwd = 2)


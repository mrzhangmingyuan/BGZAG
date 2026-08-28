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
}

# ---- 7. 输出 ----
cat("\n对称ZIP模型参数估计结果：\n")
cat(sprintf("omega_zip = %.4f (SE = %.4f)\n", result$par[1], se[1]))
cat(sprintf("alpha     = %.4f (SE = %.4f)\n", result$par[2], se[2]))
cat(sprintf("beta      = %.4f (SE = %.4f)\n", result$par[3], se[3]))
cat(sprintf("gamma     = %.4f (SE = %.4f)\n", result$par[4], se[4]))
cat(sprintf("delta     = %.4f (SE = %.4f)\n", result$par[5], se[5]))
cat(sprintf("omega     = %.4f (SE = %.4f)\n", result$par[6], se[6]))

n <- length(X_data)  
k <- length(result$par)  
log_likelihood <- -result$value  

AIC <- 2 * k - 2 * log_likelihood
BIC <- k * log(n) - 2 * log_likelihood

cat(sprintf("AIC  = %.4f\n", AIC))
cat(sprintf("BIC  = %.4f\n", BIC))
cat(sprintf("对数似然值 = %.2f\n", log_likelihood))

# ---- 8. 预测 ----
set.seed(1234)
lambda <- numeric(n)
lambda[1] <- ifelse(1 - result$par[3] != 0, result$par[6] / (1 - result$par[3]), 1e-8)
for (t in 2:n) {
  term1 <- abs(X_data[t - 1]) - result$par[4] * X_data[t - 1]
  if (term1 < 0) term1 <- 1e-8
  lambda[t] <- (result$par[6] + result$par[2] * term1^result$par[5] + 
                  result$par[3] * lambda[t - 1]^result$par[5])^(1/result$par[5])
  lambda[t] <- max(lambda[t], 1e-8)
}

predictions <- numeric(n)
lambda_prev <- lambda[1]
X_prev <- X_data[1]
for (i in 1:n) {
  term1 <- abs(X_prev) - result$par[4] * X_prev
  if (term1 < 0) term1 <- 1e-8
  lambda_pred <- (result$par[6] + result$par[2] * term1^result$par[5] + 
                    result$par[3] * lambda_prev^result$par[5])^(1/result$par[5])
  lambda_pred <- max(lambda_pred, 1e-8)
  
  if (runif(1) < result$par[1]) {
    Y_pred <- 0
  } else {
    Y_pred <- sample(c(-rpois(1, lambda_pred), rpois(1, lambda_pred)), 1)
  }
  
  predictions[i] <- Y_pred
  X_prev <- predictions[i]
  lambda_prev <- lambda_pred
}

RMS <- sqrt(mean((X_data - predictions)^2))
MAE <- mean(abs(X_data - predictions))
MSE <- mean((X_data - predictions)^2)
cat(sprintf("RMS = %.4f, MAE = %.4f, MSE = %.4f\n", RMS, MAE, MSE))

# ---- 9. PIT 和 RQR ----
lambda_hat <- lambda
pit_values <- numeric(n)

for (t in 1:n) {
  k <- X_data[t]
  
  if (k == 0) {
    zero_prob <- sym_zip_pmf(0, result$par[1], lambda_hat[t])
    pit_values[t] <- runif(1, 0, zero_prob)
  } else {
    lower <- sym_zip_cdf(k - 1, result$par[1], lambda_hat[t])
    upper <- sym_zip_cdf(k, result$par[1], lambda_hat[t])
    pit_values[t] <- runif(1, lower, upper)
  }
}

pit_values_nonzero <- pit_values[X_data != 0]

# 防止 PIT 为 0 或 1 导致 qnorm 出现 Inf
epsilon <- 1e-8
RQR <- qnorm(pmin(pmax(pit_values, epsilon), 1 - epsilon))

# ---- 绘图 ----
par(mfrow = c(1, 2))

# PIT 直方图
hist(pit_values_nonzero, breaks = 9, freq = FALSE, col = "skyblue", border = "black",
     xlab = "PIT", ylab = "Relative Frequency", main = "(c1) PIT Histogram of ZIPZAG(1,1) Model")
abline(h = 1, col = "red", lwd = 2, lty = 2)

# RQR QQ 图
qqnorm(RQR, main = "(c2) QQ Plot of RQR (ZIPZAG(1,1))")
qqline(RQR, col = "red", lwd = 2)

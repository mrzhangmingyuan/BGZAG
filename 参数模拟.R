library(stats)
library(ggplot2)
library(tidyr)
library(dplyr)
set.seed(1)
# DND PMF
dnd_pmf <- function(k, lambda) {
  pnorm((k + 1) / lambda) - pnorm(k / lambda)
}

# 采样一个 DND(0, lambda)
sample_dnd <- function(lambda, k_max = 60) {
  k_vals <- -k_max:k_max
  probs <- sapply(k_vals, dnd_pmf, lambda = lambda)
  
  probs <- probs / sum(probs)
  sample(k_vals, size = 1, prob = probs)
}
# 模拟 DND-GARCH(1,1) 时间序列
simulate_dnd_series <- function(n = 200, omega = 0.5, alpha = 0.2, beta = 0.3, gamma = -0.3) {
  X <- numeric(n)
  lambda <- numeric(n)
  lambda[1] <- ifelse(1 - beta != 0, omega / (1 - beta), 1e-8)
  
  for (t in 2:n) {
    lambda[t] <- omega + alpha * (abs(X[t - 1]) - gamma * X[t - 1])^2 + beta * lambda[t - 1]
    lambda[t] <- max(lambda[t], 1e-8)
    X[t] <- sample_dnd(lambda[t])
  }
  
  return(X)
}

# 对数似然函数
negative_log_likelihood <- function(params, X_data) {
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  gamma <- params[4]
  
  n <- length(X_data)
  lambda <- numeric(n)
  lambda[1] <- ifelse(1 - beta != 0, omega / (1 - beta), 1e-8)
  log_likelihood <- 0
  
  for (t in 2:n) {
    lambda[t] <- omega + alpha * (abs(X_data[t - 1]) - gamma * X_data[t - 1])^2 + beta * lambda[t - 1]
    lambda[t] <- max(lambda[t], 1e-8)
    
    k_vals <- -60:60
    pmf_vals <- sapply(k_vals, dnd_pmf, lambda = lambda[t])
    pmf_vals <- pmf_vals / sum(pmf_vals)
    
    idx <- which(k_vals == X_data[t])
    if (length(idx) == 0) next
    prob <- max(pmf_vals[idx], 1e-8)
    log_likelihood <- log_likelihood + log(prob)
  }
  
  return(-log_likelihood)
}

# 参数初始化
params_init <- c(0.5, 0.2, 0.3, -0.3)
true_params <- c(omega = 0.5, alpha = 0.2, beta = 0.3, gamma = -0.3)

# 不同样本量
sample_sizes <- c(500, 1000, 1500)
num_iterations <- 500  # 可设置为500

# 存储所有结果
all_results <- list()

# 外层循环：遍历不同样本量
for (n in sample_sizes) {
  cat(sprintf("\n============ 开始样本量 n = %d 的模拟 ============\n", n))
  
  estimated_params_list <- list()
  
  for (i in 1:num_iterations) {
    X_data <- simulate_dnd_series(n = n)
    
    res <- tryCatch({
      optim(
        par = params_init,
        fn = negative_log_likelihood,
        X_data = X_data,
        method = "L-BFGS-B",
        lower = c(1e-8, 0, 0, -1),
        upper = c(5, 1, 1, 1)
      )
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$convergence == 0) {
      estimated_params_list[[length(estimated_params_list) + 1]] <- res$par
    } else {
      cat(sprintf("第 %d 次优化失败,跳过。\n", i))
    }
  }
  
  # 统计结果
  if (length(estimated_params_list) > 0) {
    estimated_params_matrix <- do.call(rbind, estimated_params_list)
    mean_params <- colMeans(estimated_params_matrix)
    SEs <- apply(estimated_params_matrix, 2, sd) / sqrt(nrow(estimated_params_matrix))
    MADEs <- colMeans(abs(sweep(estimated_params_matrix, 2, mean_params)))
    
    # 输出结果
    cat(sprintf("\n完成 %d 次模拟与估计 (n = %d)\n", length(estimated_params_list), n))
    cat("平均估计参数：\n")
    cat(sprintf("omega: %.4f    SE: %.4f    MADE: %.4f\n", mean_params[1], SEs[1], MADEs[1]))
    cat(sprintf("alpha: %.4f    SE: %.4f    MADE: %.4f\n", mean_params[2], SEs[2], MADEs[2]))
    cat(sprintf("beta : %.4f    SE: %.4f    MADE: %.4f\n", mean_params[3], SEs[3], MADEs[3]))
    cat(sprintf("gamma: %.4f    SE: %.4f    MADE: %.4f\n", mean_params[4], SEs[4], MADEs[4]))
    
    # 存储结果用于箱线图
    df <- data.frame(estimated_params_matrix)
    colnames(df) <- c("omega", "alpha", "beta", "gamma")
    df$sample_size <- as.factor(n)
    all_results[[length(all_results) + 1]] <- df
  }
}

# 合并所有结果
combined_results <- do.call(rbind, all_results)

# 转换为长格式用于绘图
long_results <- combined_results %>%
  pivot_longer(cols = c("omega", "alpha", "beta", "gamma"),
               names_to = "parameter",
               values_to = "estimate")

# 创建箱线图
true_value_df <- data.frame(
  parameter = c("omega", "alpha", "beta", "gamma"),
  true_value = true_params
)

p <- ggplot(long_results, aes(x = sample_size, y = estimate, fill = sample_size)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("#317CB7", "#6DADDE", "#B6D7E8")) +
  geom_hline(data = true_value_df, 
             aes(yintercept = true_value), 
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_text(data = true_value_df,
            aes(x = Inf, y = true_value, 
                label = paste0("True = ", true_value)),
            hjust = 1.05, vjust = -0.5, 
            color = "red", size = 3.5, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~parameter, scales = "free_y", ncol = 2) +
  labs(
    title = "(a) Box Plot of Parameter Set A1",
    x = "n",
    y = "Estimated Value",
    fill = "n"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom"
  )

print(p)

# 输出汇总统计表
cat("\n\n============ 汇总统计表 ============\n")
summary_stats <- long_results %>%
  group_by(sample_size, parameter) %>%
  summarise(
    Mean = mean(estimate),
    SD = sd(estimate),
    Bias = mean(estimate) - true_params[parameter[1]],
    .groups = "drop"
  )
print(summary_stats)
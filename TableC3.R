# ==============================================================================
# MSc Final Project
# Reliability Analysis of Carbon-Film Resistors (Two-Stage Method)
# Reference: Meeker & Escobar (1998)
# Data Source: SMRD.data::resistor (Table C.3)
# ==============================================================================

library(tidyverse)
library(survival)
library(SMRD.data)
library(ggfortify)

# ==============================================================================
# Step 1: Data Loading and Preprocessing
# ==============================================================================

data("resistor")

df <- resistor %>%
  rename(
    Unit = resistor,
    Temp = celsius,
    Time = hours,
    Pct_Increase = percent
  )

# Correct unit indexing (Unit 10 does not exist)
correct_units <- c(1:9, 11:30)
df$Unit <- rep(correct_units, each = 4)
df$Temp <- c(rep(83, 9*4), rep(133, 10*4), rep(173, 10*4))

# ==============================================================================
# Step 2: Degradation Modeling (Stage 1)
# Soft Failure Threshold: 5%
# ==============================================================================

Df <- 5  # 5% resistance increase defines failure

# Spaghetti Plot
ggplot(df, aes(Time, Pct_Increase, group = Unit, color = factor(Temp))) +
  geom_line() + geom_point() +
  geom_hline(yintercept = Df, linetype = "dashed", color = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Degradation Paths of Carbon-Film Resistors",
       x = "Time (Hours, log scale)",
       y = "Resistance Increase (%)",
       color = "Temperature (°C)") +
  theme_bw()

# ---- Pseudo-Failure Time Estimation ----
# Power Law: log(D) = beta0 + beta1 * log(t)
# Only beta1 > 0 is physically meaningful

pseudo_failures <- df %>%
  filter(Pct_Increase > 0) %>%   # avoid log(0)
  group_by(Unit, Temp) %>%
  do({
    m <- try(lm(log(Pct_Increase) ~ log(Time), data = .), silent = TRUE)
    if (inherits(m, "try-error")) return(data.frame(T_fail = NA))
    
    beta0 <- coef(m)[1]
    beta1 <- coef(m)[2]
    
    if (beta1 <= 0) return(data.frame(T_fail = NA))  # non-degrading unit
    
    T_fail <- exp((log(Df) - beta0) / beta1)
    data.frame(T_fail = T_fail)
  }) %>%
  ungroup() %>%
  filter(!is.na(T_fail)) %>%
  mutate(
    TempK = Temp + 273.15,
    Arrhenius_x = 11605 / TempK
  )

print(pseudo_failures)

# ==============================================================================
# Step 3: Weibull–Arrhenius Acceleration Model (Stage 2)
# survreg fits SEV on log(T)
# ==============================================================================

model_mle <- survreg(
  Surv(T_fail) ~ Arrhenius_x,
  data = pseudo_failures,
  dist = "weibull"
)

summary(model_mle)

mu_hat    <- coef(model_mle)[1]
beta_hat  <- coef(model_mle)[2]   # proportional to Ea / k
sigma_hat <- model_mle$scale      # Weibull scale (1 / shape)

# NOTE:
# beta_hat is proportional to activation energy (Ea/k), not Ea itself

# ==============================================================================
# Step 4: Prediction at Use Condition (50°C)
# ==============================================================================

use_temp_C <- 50
use_x <- 11605 / (use_temp_C + 273.15)

mu_use <- mu_hat + beta_hat * use_x

# Mean Time To Failure
MTTF <- exp(mu_use) * gamma(1 + sigma_hat)

# B10 Life (p = 0.1)
p <- 0.1
z_p <- log(-log(1 - p))
B10 <- exp(mu_use + sigma_hat * z_p)

# Failure probability at 100,000 hours
t_eval <- 1e5
z_eval <- (log(t_eval) - mu_use) / sigma_hat
F_t <- 1 - exp(-exp(z_eval))

cat("--- Predictions at 50°C ---")
cat("MTTF:", MTTF)
cat("B10 Life:", B10)
cat("Failure Prob at 100k hours:", F_t)

# ==============================================================================
# Step 5: Normal Approximation CI for B10
# (Asymptotic, conditional on sigma)
# ==============================================================================

pred <- predict(
  model_mle,
  newdata = data.frame(Arrhenius_x = use_x),
  type = "lp",
  se.fit = TRUE
)

# B10 quantile
p <- 0.1
z_p <- log(-log(1 - p))

# Linear predictor
eta_hat <- pred$fit

# log B10 estimate
log_B10_hat <- eta_hat + sigma_hat * z_p

# Variance–covariance matrix
V <- vcov(model_mle)

# Gradient
grad <- c(
  1,
  use_x,
  sigma_hat * z_p
)

# Delta-method SE
se_log_B10 <- sqrt(t(grad) %*% V %*% grad)

# 95% Normal CI
log_B10_lower <- log_B10_hat - 1.96 * se_log_B10
log_B10_upper <- log_B10_hat + 1.96 * se_log_B10

# Back-transform
B10_hat   <- exp(log_B10_hat)
B10_lower <- exp(log_B10_lower)
B10_upper <- exp(log_B10_upper)

cat("B10 =", B10_hat)
cat("95% Normal CI:", B10_lower, "-", B10_upper, "\n")


# ==============================================================================
# Step 6: Bootstrap (Two-Stage Uncertainty)
# ==============================================================================
library(boot)

# 1. 定義 Bootstrap 估計函數 (這是為了符合 boot 套件格式)
boot_func <- function(data, indices) {
  # 重抽樣 Units
  units <- unique(data$Unit)
  resampled_units <- units[indices]
  
  # 建立對應的數據集
  boot_df <- map_dfr(seq_along(resampled_units), function(j) {
    data %>% filter(Unit == resampled_units[j]) %>% mutate(NewID = j)
  })
  
  # 第一階段：估計虛擬失效時間
  boot_fail <- boot_df %>%
    group_by(NewID, Temp) %>%
    do({
      m <- try(lm(log(Pct_Increase) ~ log(Time), data = .), silent = TRUE)
      if (inherits(m, "try-error")) return(data.frame(T_fail = NA))
      b1 <- coef(m)[2]
      if (b1 <= 0) return(data.frame(T_fail = NA))
      data.frame(T_fail = exp((log(Df) - coef(m)[1]) / b1))
    }) %>%
    ungroup() %>% filter(!is.na(T_fail)) %>%
    mutate(Arrhenius_x = 11605 / (Temp + 273.15))
  
  # 第二階段：擬合加速模型
  fit <- try(survreg(Surv(T_fail) ~ Arrhenius_x, data = boot_fail, dist = "weibull"), silent = TRUE)
  
  if (inherits(fit, "try-error")) return(NA)
  
  # 計算 B10 (務必使用正確的 z_p!)
  z_p_correct <- log(-log(1 - 0.1))
  b10_est <- exp(coef(fit)[1] + coef(fit)[2] * use_x + fit$scale * z_p_correct)
  return(log(b10_est)) # 返回對數值，分佈會較對稱
}

# 2. 執行 boot
# 注意：為了計算 BCa，R 的 boot 函數對單純重抽樣的支援較複雜，
# 這裡我們使用較簡單的 percentile 呈現，並附上偏差修正的觀念
set.seed(5014)
# 由於 BCa 需要 Jackknife (逐一剔除法) 來計算加速度 a，
# 在兩階段模型中運算量極大且容易不收斂。
# 建議簡報中展示 "Bias-Corrected" 的概念：

# 計算偏差 (Bias)
mle_log_b10 <- log(B10) # 原始 MLE 的對數 B10
z0 <- qnorm(sum(log(boot_results) < mle_log_b10) / length(boot_results))
cat("Bias correction factor z0:", z0)

# 如果 z0 遠離 0，代表有顯著偏差。

# 假設你的 boot 執行結果儲存在 res_boot 中
res_boot <- boot(data = df, statistic = boot_func, R = 1000)

# 1. 獲取 BCa 區間 (這會自動包含 z0 修正)
boot_ci_all <- boot.ci(res_boot, type = c("perc", "bca"))
print(boot_ci_all)

# 2. 轉換回原始時間尺度 (因為 boot_func 回傳的是 log)
bca_lower <- exp(boot_ci_all$bca[4])
bca_upper <- exp(boot_ci_all$bca[5])
perc_lower <- exp(boot_ci_all$percent[4])
perc_upper <- exp(boot_ci_all$percent[5])

# 3. 點估計 (使用偏差修正後的點估計)
# 在 Bootstrap 中，修正後的點估計通常透過 bias 調整
b10_log_mean <- mean(res_boot$t)
b10_point_corrected <- exp(mle_log_b10 - (b10_log_mean - mle_log_b10))

cat("--- Bootstrap 最終結果 (B10) ---")
cat("修正後點估計:", b10_point_corrected)
exp(median(res_boot$t))
cat("BCa 95% 區間:", bca_lower, "-", bca_upper)
cat("Percentile 95% 區間:", perc_lower, "-", perc_upper)

# ==============================================================================
# Step 7: Reverse Prediction – Max Allowable Temperature
# ==============================================================================

target_life <- 43800  # 5 years
numerator <- log(target_life) - mu_hat - sigma_hat * z_p
arr_x_needed <- numerator / beta_hat

max_temp_C <- 11605 / arr_x_needed - 273.15

cat("Max allowable temperature for B10 ≥ 5 years:",
    round(max_temp_C, 2), "°C\n")

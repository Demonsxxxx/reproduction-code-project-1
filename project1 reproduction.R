## ---- RQ5: one-year-ahead forecast with metro-specific accuracy ----
set.seed(123)

# 1) Metro-specific accuracy means (approx from paper's Total/HUD ratios)
mu_map <- c(
  "New York"         = 0.926,  # 73,523 / 79,348
  "Los Angeles"      = 0.715,  # 46,874 / 65,585
  "Washington, D.C." = 0.958,  # 8,350 / 8,719
  "Seattle"          = 0.824,  # 10,730 / 13,027
  "Houston"          = 0.724   # 4,031 / 5,565
)
var_pi   <- 0.0015         # prior variance used in the paper examples
sigma2_psi <- 0.001        # small innovation variance for log-odds
phi_bar <- 0.94            # pooled rent effect prior mean (paper Sec. 4.2)
DeltaZRI_default <- 0.05   # assume +5% YoY ZRI (edit with real values if you have them)

# HUD counts (last observed year)
hud_counts <- c(
  "New York"=73523,
  "Los Angeles"=46874,
  "Washington, D.C."=8350,
  "Seattle"=10730,
  "Houston"=4031
)

# ---- helpers ----
beta_from_mean_var <- function(mu, v){
  k <- mu*(1-mu)/v - 1
  c(a = mu*k, b = (1-mu)*k)
}
log_beta_binom <- function(C, H, a, b){
  lchoose(H, C) + lbeta(C + a, H - C + b) - lbeta(a, b)
}
posterior_H_given_C <- function(C, mu, v=var_pi, H_mult=1.6){
  ab <- beta_from_mean_var(mu, v); a <- ab["a"]; b <- ab["b"]
  H_grid <- C:ceiling(H_mult*C)
  loglik <- vapply(H_grid, function(H) log_beta_binom(C, H, a, b), numeric(1))
  w <- exp(loglik - max(loglik)); w <- w / sum(w)
  list(H_grid=H_grid, w=w)
}
forecast_one_step <- function(metro, C, S=6000, dZRI=DeltaZRI_default){
  mu <- mu_map[[metro]]
  post <- posterior_H_given_C(C, mu, var_pi)
  H_t <- sample(post$H_grid, size=S, replace=TRUE, prob=post$w)
  
  # working denominator; if you know N_t, plug it in here instead
  N_t <- pmax(round(H_t / mu), H_t + 1L)
  p_t <- pmin(pmax(H_t / N_t, 1e-6), 1 - 1e-6)
  psi_t <- qlogis(p_t)
  
  # evolve one step with positive rent shock
  psi_tp1 <- psi_t + phi_bar * dZRI + rnorm(S, 0, sqrt(sigma2_psi))
  p_tp1 <- plogis(psi_tp1)
  H_tp1 <- rbinom(S, size=N_t, prob=p_tp1)
  
  mean_true <- mean(H_tp1)
  lo_true   <- unname(quantile(H_tp1, 0.025))
  hi_true   <- unname(quantile(H_tp1, 0.975))
  prob_inc  <- mean(H_tp1 > H_t)
  
  # also show expected observed (HUD-like) counts using next-year E[pi]
  mean_obs <- mean_true * mu
  lo_obs   <- lo_true * mu
  hi_obs   <- hi_true * mu
  
  c(mean_true=mean_true, lo_true=lo_true, hi_true=hi_true,
    mean_obs=mean_obs, lo_obs=lo_obs, hi_obs=hi_obs,
    prob_inc=prob_inc)
}

# run for all metros
res <- t(mapply(forecast_one_step, names(hud_counts), hud_counts))
tab <- data.frame(
  Metro = names(hud_counts),
  `HUD count (C_t)`     = as.integer(hud_counts),
  `Forecast TRUE mean`  = as.integer(round(res[,"mean_true"])),
  `TRUE 95% lo`         = as.integer(round(res[,"lo_true"])),
  `TRUE 95% hi`         = as.integer(round(res[,"hi_true"])),
  `Forecast OBS mean`   = as.integer(round(res[,"mean_obs"])),
  `Prob(increase)`      = round(res[,"prob_inc"], 3),
  check.names = FALSE
)
print(tab, row.names = FALSE)

# Tip: if you have metro-specific ΔZRI, pass a named vector and call:
# res <- t(mapply(forecast_one_step, names(hud_counts), hud_counts, dZRI = DeltaZRI_by_metro))

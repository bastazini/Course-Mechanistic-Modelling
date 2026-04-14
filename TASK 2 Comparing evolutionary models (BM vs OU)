################################################################################
##                                                                            ##
##  TASK 2                                                                    ##
##  Brownian Motion vs Ornstein-Uhlenbeck -- comparing evolutionary models   ##
##                                                                            ##
##  BACKGROUND                                                                ##
##  BM: trait increments are i.i.d. N(0, sigma2). Variance grows without     ##
##  bound as sigma2 * t. The trajectory is a random walk with no attractor.  ##
##                                                                            ##
##  OU: a deterministic pull towards optimum theta is added:                 ##
##    X(t) = X(t-1) + alpha*(theta - X(t-1)) + N(0, sigma2)                 ##
##  alpha controls the strength of stabilising selection. Large alpha ->      ##
##  strong pull, trait stays near theta. alpha = 0 recovers BM.              ##
##                                                                            ##
##  We simulate 10 BM trajectories with increasing sigma2 (red palette) and   ##
##  10 OU trajectories with increasing alpha at fixed sigma2 (blue palette).  ##
##                                                                            ##
################################################################################


# -------------------------------
# 0. Clean workspace
# -------------------------------
rm(list = ls())


# -------------------------------
# 1. Set random seed (reproducibility)
# -------------------------------
set.seed(12)


# -------------------------------
# 2. Define simulation parameters
# -------------------------------

n_steps <- 200
n_traj  <- 200       # trajectories per model
theta   <- 0         # OU optimum


# -------------------------------
# 3. Define parameter gradients
# -------------------------------

# BM: increasing variance per step
sigma2_bm <- seq(0.02, 0.15, length.out = n_traj)

# OU: increasing selection strength
alpha_ou  <- seq(0.05, 0.3, length.out = n_traj)

# OU: fixed noise variance
sigma2_ou <- 0.1


# -------------------------------
# 4. Initialize storage matrix
# -------------------------------

# Columns 1:n_traj = BM
# Columns (n_traj+1):(2*n_traj) = OU
traj <- matrix(0, nrow = n_steps, ncol = n_traj * 2)


# -------------------------------
# 5. Simulate BM and OU processes
# -------------------------------

for (i in seq_len(n_traj)) {
  
  # ---- BM process ----
  # Cumulative sum of N(0, sigma2) increments
  traj[, i] <- cumsum(rnorm(n_steps, mean = 0, sd = sqrt(sigma2_bm[i])))
  
  
  # ---- OU process ----
  # Euler-Maruyama discretisation, starting at 0
  x <- 0  
  traj[1, i + n_traj] <- x
  
  for (t in 2:n_steps) {
    # Deterministic pull + stochastic noise
    x <- x + alpha_ou[i] * (theta - x) +
      rnorm(1, mean = 0, sd = sqrt(sigma2_ou))
    
    traj[t, i + n_traj] <- x
  }
}


# -------------------------------
# 6. Define color palettes
# -------------------------------

# Light = weak parameter, Dark = strong parameter
cols_bm <- colorRampPalette(c("#fcbba1", "#a50f15"))(n_traj)  # BM (red)
cols_ou <- colorRampPalette(c("#c6dbef", "#08306b"))(n_traj)  # OU (blue)


# -------------------------------
# 7. Plot trajectories
# -------------------------------

matplot(traj,
        type = "l",
        col  = c(cols_bm, cols_ou),
        lwd  = 2,
        lty  = 1,
        xlab = "Time (generations)",
        ylab = "Trait value",
        main = "Task 2 -- BM (red) vs OU (blue): trait trajectories")


# -------------------------------
# 8. Add OU optimum reference line
# -------------------------------

abline(h = theta, lty = 2, col = "black", lwd = 1)


# -------------------------------
# 9. Label the optimum
# -------------------------------

text(x = n_steps * 0.98, y = theta + 0.05,
     labels = paste0("optimum theta = ", theta),
     adj = 1, cex = 0.8, col = "grey30")


# -------------------------------
# 10. Add legend
# -------------------------------

legend("topleft",
       legend = c(paste0("BM  sigma2 = ", round(sigma2_bm, 3)),
                  paste0("OU  alpha = ",  round(alpha_ou,  2))),
       col    = c(cols_bm, cols_ou),
       lty    = 1, lwd = 2, bty = "n", cex = 0.75,
       title  = "Model and parameter")

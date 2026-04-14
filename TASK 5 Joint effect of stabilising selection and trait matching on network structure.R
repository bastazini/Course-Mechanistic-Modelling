################################################################################
##                                                                            ##
##  TASK 5                                                                    ##
##  Joint effect of stabilising selection and trait matching on network       ##
##  structure                                                                 ##
##                                                                            ##
##  BACKGROUND                                                               ##
##  Two parameters jointly determine the structure of an STC network:        ##
##                                                                           ##
##  1. STABILISING SELECTION (alpha, the OU pull parameter)                  ##
##     Under the OU model at stationarity, trait values are drawn from:      ##
##       X ~ Normal(theta,  sigma2 / (2 * alpha))                            ##
##     High alpha -> small variance -> all species have similar traits ->    ##
##     many trait pairs overlap -> dense, non-selective network.             ##
##     Low alpha -> large variance -> traits spread widely -> only closely   ##
##     matched species interact -> sparse, modular network.                  ##
##     Note: alpha = 0 (pure BM) has no stationary distribution and is       ##
##     therefore excluded from this task.                                    ##
##                                                                           ##
##  2. TRAIT MATCHING (delta_max, upper bound of per-species trait range)    ##
##     In the STC rule, species i and j interact if:                         ##
##       |V_i - W_j|  <  0.5 * (delta_L_i + delta_H_j)                       ##
##     Each species is assigned a range delta drawn from                     ##
##       Uniform(delta_min, delta_max)                                       ##
##     Wide delta -> generalist (interacts with many partners).              ##
##     Narrow delta -> specialist (interacts only with close trait matches). ##
##     delta_max sets the ceiling of this variability across the community.  ##
##                                                                           ##
##  AVAILABLE METRICS (set metric_name below)                                ##
##  Any metric accepted by bipartite::networklevel(index = ...):             ##
##    "connectance"        -- proportion of realised links                   ##
##    "nestedness"         -- NODF nestedness                                ##
##    "modularity"         -- Q modularity                                   ##
##    "H2"                 -- network-level specialisation (0 = generalist)  ##
##    "links per species"  -- mean degree across both guilds                 ##
##  Run networklevel(web, index = "ALLBUTDD") on a test web to see all.      ##
##                                                                           ##
##                                                                            ##
################################################################################

# -------------------------------
# 0. Clean workspace
# -------------------------------
rm(list = ls())


# -------------------------------
# 1. Load required packages
# -------------------------------
library(ape)        
library(bipartite)  
library(ggplot2)    

# -------------------------------
# 2. Set seed
# -------------------------------
set.seed(123)

# -------------------------------
# 3. Parameters
# -------------------------------

metric_name <- "modularity"
# Network-level property to evaluate.
# "links per species" captures overall interaction density.

n_L <- 20  # consumer species (rows)
n_H <- 20  # resource (columns)
# Defines bipartite network size. Larger values increase realism but cost.

n_rep <- 100
# Number of stochastic replicates per parameter combination.
# Higher values → smoother heatmaps, lower variance.

sigma2 <- 1.0
# Brownian motion variance component.
# Controls evolutionary noise intensity in OU process.

theta_L <- 0
theta_H <- 0
# OU optima: both guilds evolve toward same adaptive peak.
# Setting them equal maximises trait overlap potential.

alpha_grid <- c(0.1, 0.5, 1, 2, 4, 8)
# Strength of stabilising selection:
# low alpha → wide trait dispersion
# high alpha → strong constraint around theta (trait clustering)

delta_max_grid <- c(0.1, 0.3, 0.6, 0.9, 1.2)
delta_min <- 0.05
# STC trait tolerance:
# each species has its own interaction “width”
# small values → specialists
# large values → generalists

# -------------------------------
# 4. Phylogenies (fixed across simulations)
# -------------------------------
tree_L <- rcoal(n_L)
tree_H <- rcoal(n_H)

tree_L$tip.label <- sprintf("L_%03d", seq_len(n_L))
tree_H$tip.label <- sprintf("H_%03d", seq_len(n_H))

# Coalescent trees simulate shared evolutionary history.
# Fixing trees across runs isolates effects of alpha and delta.

# -------------------------------
# 5. Simulation function
# -------------------------------
simulate_metric <- function(alpha_L, alpha_H,
                            delta_min, delta_max,
                            metric) {
  
  # -------------------------------
  # 5.1 Trait evolution (OU model on phylogeny)
  # -------------------------------
  trait_L <- rTraitCont(tree_L,
                        model = if (alpha_L == 0) "BM" else "OU",
                        alpha = alpha_L,
                        sigma = sqrt(sigma2),
                        theta = theta_L)
  
  trait_H <- rTraitCont(tree_H,
                        model = if (alpha_H == 0) "BM" else "OU",
                        alpha = alpha_H,
                        sigma = sqrt(sigma2),
                        theta = theta_H)
  
  # Each species inherits a trait value evolving along its phylogeny.
  # OU introduces stabilising selection:
  #   - high alpha → convergence to theta
  #   - low alpha → high variance among species
  
  # -------------------------------
  # 5.2 Species-specific trait tolerance (STC variability)
  # -------------------------------
  delta_L <- runif(n_L, delta_min, delta_max)
  delta_H <- runif(n_H, delta_min, delta_max)
  
  # Each species has a trait "interaction width".
  # This introduces ecological heterogeneity:
  #   - small delta → specialist
  #   - large delta → generalist
  
  # -------------------------------
  # 5.3 Interaction matrix (Santamaría STC rule)
  # -------------------------------
  web <- outer(seq_len(n_L), seq_len(n_H),
               Vectorize(function(i, j)
                 as.integer(
                   abs(trait_L[i] - trait_H[j]) <=
                     0.5 * (delta_L[i] + delta_H[j])
                 )
               ))
  
  # Core ecological assumption:
  # interaction occurs if trait distributions overlap.
  #
  # This is equivalent to:
  # |mean difference| <= half-sum of trait ranges
  #
  # → species interact if their “niches overlap”
  
  # -------------------------------
  # 5.4 Ensure no isolated species
  # -------------------------------
  for (j in which(colSums(web) == 0))
    web[which.min(abs(trait_L - trait_H[j])), j] <- 1L
  
  for (i in which(rowSums(web) == 0))
    web[i, which.min(abs(trait_H - trait_L[i]))] <- 1L
  
  # Ecological constraint:
  # every species must have at least one interaction
  # → avoids degenerate rows/columns in network metrics
  
  # -------------------------------
  # 5.5 Compute network metric
  # -------------------------------
  val <- tryCatch(
    networklevel(web, index = metric),
    error = function(e) NA_real_
  )
  
  # Defensive programming:
  # prevents simulation failure from stopping full grid run
  
  as.numeric(val[1])
}

# -------------------------------
# 6. Grid simulation
# -------------------------------
cat(sprintf("Running %d simulations...\n",
            length(alpha_grid) * length(delta_max_grid) * n_rep))

grid_results <- do.call(rbind, lapply(alpha_grid, function(a) {
  do.call(rbind, lapply(delta_max_grid, function(dmax) {
    
    vals <- replicate(n_rep,
                      simulate_metric(
                        alpha_L = a,
                        alpha_H = a,
                        delta_min = delta_min,
                        delta_max = dmax,
                        metric = metric_name
                      ))
    
    data.frame(
      alpha     = a,
      delta_max = dmax,
      mean_val  = mean(vals, na.rm = TRUE),
      se_val    = sd(vals, na.rm = TRUE) / sqrt(sum(!is.na(vals)))
    )
    
    # Aggregation:
    # mean → expected network structure
    # SE   → uncertainty across stochastic replicates
  }))
}))

cat("Done.\n")

# -------------------------------
# 7. Format axes for plotting
# -------------------------------
grid_results$alpha_label <- factor(grid_results$alpha,
                                   levels = alpha_grid)

grid_results$delta_label <- factor(grid_results$delta_max,
                                   levels = delta_max_grid)

# Ensures correct ordering in heatmap (not alphabetical)

# -------------------------------
# 8. Heatmap visualization
# -------------------------------
print(
  ggplot(grid_results,
         aes(x = delta_label,
             y = alpha_label,
             fill = mean_val)) +
    
    geom_tile(color = "white") +
    # Grid cells represent parameter combinations
    
    geom_text(aes(label = round(mean_val, 2)),
              size = 4, color = "white") +
    # Overlay numeric values for interpretability
    
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = median(grid_results$mean_val, na.rm = TRUE),
      name = metric_name
    ) +
    # Diverging scale highlights low vs high connectivity regimes
    
    labs(
      title = "TASK 5 - Joint effect of stabilising selection and trait matching",
      x = "Trait tolerance (delta_max)",
      y = "Stabilising selection strength (alpha)"
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
)

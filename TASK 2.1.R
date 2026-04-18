# ================================
#TASK 2.1 
#BM vs OU trait evolution on same tree
# Base R + ape + phytools only
# ================================
rm(list = ls())

library(ape)
library(phytools)

set.seed(123)

# -------------------------------
# 1. Simulate tree
# -------------------------------
tree <- rcoal(40)

# -------------------------------
# 2. Simulate traits using rTraitCont
# -------------------------------
traits_BM <- rTraitCont(tree, model = "BM", sigma2 = 0.1)

traits_OU <- rTraitCont(tree, model = "OU",
                        alpha = 10,
                        sigma2 = 0.1,
                        theta = 0)

# -------------------------------
# 3 Ploting
# -------------------------------

# ===============================
# PANEL 1: Density comparison
# ===============================
dens_BM <- density(traits_BM)
dens_OU <- density(traits_OU)

plot(dens_BM,
     col = "blue",
     lwd = 2,
     main = "Trait distributions",
     xlab = "Trait value",
     ylim = range(c(dens_BM$y, dens_OU$y)))

lines(dens_OU, col = "red", lwd = 2)

legend("topright",
       legend = c("BM", "OU"),
       col = c("blue", "red"),
       lwd = 2,
       bty = "n")


### set pannels
par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))

# ===============================
# PANEL 2: Mapping BM  continuous character.
# ===============================
obj_BM <- contMap(tree, traits_BM, plot = FALSE)
plot(obj_BM, main = "BM on phylogeny")

# ===============================
# PANEL 3: Mapping BM  continuous character.
# ===============================
obj_OU <- contMap(tree, traits_OU, plot = FALSE)
plot(obj_OU, main = "OU on phylogeny")

# ===============================
# PANEL 4: Tip-state comparison
# ===============================

plot(tree,
     show.tip.label = FALSE,
     main = "Tip states (BM)")

tiplabels(pch = 16,
          col = heat.colors(50)[rank(traits_BM)])


plot(tree,
     show.tip.label = FALSE,
     main = "Tip states (OU)")

tiplabels(pch = 16,
          col = heat.colors(50)[rank(traits_OU)])


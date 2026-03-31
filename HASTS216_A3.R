# ============================================================================
# HASTS 416 Tutorial – Complete R Solution
# ============================================================================

library(markovchain)
library(igraph)
library(ggplot2)
library(reshape2)




# Question 3 – Traffic model
# ============================================================================
cat("\n\n")
cat("==========================================\n")
cat("Question 3: Traffic Model\n")
cat("==========================================\n\n")

# Define state space
states <- c("light", "heavy", "jammed")

# Transition matrices
P1 <- matrix(c(0.4, 0.4, 0.2,
               0.3, 0.4, 0.2,
               0.0, 0.1, 0.9), nrow = 3, byrow = TRUE)

P2 <- matrix(c(0.1, 0.5, 0.4,
               0.1, 0.3, 0.6,
               0.0, 0.1, 0.9), nrow = 3, byrow = TRUE)

rownames(P1) <- rownames(P2) <- states
colnames(P1) <- colnames(P2) <- states

cat("Period 1 (1PM-4PM) Transition Matrix:\n")
print(P1)
cat("\nPeriod 2 (4PM-6PM) Transition Matrix:\n")
print(P2)

#> cat("Period 1 (1PM-4PM) Transition Matrix:\n")
#Period 1 (1PM-4PM) Transition Matrix:
#  > print(P1)
#light heavy jammed
#light    0.4   0.4    0.2
#heavy    0.3   0.4    0.2
#jammed   0.0   0.1    0.9
#> cat("\nPeriod 2 (4PM-6PM) Transition Matrix:\n")

#Period 2 (4PM-6PM) Transition Matrix:
#  > print(P2)
#light heavy jammed
#light    0.1   0.5    0.4
#heavy    0.1   0.3    0.6
#jammed   0.0   0.1    0.9


steps1 <- 9
steps2 <- 6

cat("\nTime steps:\n")
cat("  • 1PM to 4PM:", steps1, "steps (20 min each)\n")
cat("  • 4PM to 6PM:", steps2, "steps (20 min each)\n")



#Time steps:
#  > cat("  • 1PM to 4PM:", steps1, "steps (20 min each)\n")
#• 1PM to 4PM: 9 steps (20 min each)
#> cat("  • 4PM to 6PM:", steps2, "steps (20 min each)\n")
#• 4PM to 6PM: 6 steps (20 min each)
 


# ============================================================================
# (a) Analytical distribution
# ============================================================================
cat("\n--- (a) Analytical Distribution ---\n")

dist <- c(1, 0, 0)
names(dist) <- states
cat("Initial distribution (1PM):", round(dist, 4), "\n")

for(i in 1:steps1) {
  dist <- dist %*% P1
}
cat("\nDistribution at 4PM:", round(dist, 4), "\n")

for(i in 1:steps2) {
  dist <- dist %*% P2
}
cat("\nDistribution at 6PM:", round(dist, 4), "\n")

# Distribution at 6PM: 0.0118 0.1059 0.6805 



# ============================================================================
# (b) Simulation
# ============================================================================
cat("\n--- (b) Simulation Verification ---\n")

set.seed(456)
n_sim <- 10000
final_states <- character(n_sim)

for(i in 1:n_sim) {
  state <- "light"
  
  for(step in 1:steps1) {
    probs <- P1[state, ]
    state <- sample(states, size = 1, prob = probs)
  }
  
  for(step in 1:steps2) {
    probs <- P2[state, ]
    state <- sample(states, size = 1, prob = probs)
  }
  
  final_states[i] <- state
}

sim_dist <- table(final_states) / n_sim
cat("\nSimulated distribution (10,000 trajectories):\n")
print(round(sim_dist, 4))

#> print(round(sim_dist, 4))
#final_states
#heavy jammed  light 
#0.1326 0.8535 0.0139 


# Convert sim_dist to a vector in the same order as dist
sim_dist_vector <- as.vector(sim_dist[states])
names(sim_dist_vector) <- states

# Comparison
cat("\nComparison:\n")
comparison <- rbind(Analytical = dist, Simulated = sim_dist_vector)
print(round(comparison, 4))

#Comparison:
#  > comparison <- rbind(Analytical = dist, Simulated = sim_dist_vector)
#> print(round(comparison, 4))
#light  heavy jammed
#0.0118 0.1059 0.6805
#Simulated 0.0139 0.1326 0.8535


# Difference - FIXED: use aligned vectors
diff <- abs(dist - sim_dist_vector)
cat("\nAbsolute difference:\n")
print(round(diff, 4))


#Absolute difference:
#  > print(round(diff, 4))
#light  heavy jammed
#[1,] 0.0021 0.0267  0.173



# Chi-square test
expected <- dist * n_sim
observed <- as.vector(table(factor(final_states, levels = states)))
chi_sq <- sum((observed - expected)^2 / expected)
p_value <- 1 - pchisq(chi_sq, df = 2)



cat("\nChi-square goodness-of-fit test:\n")
cat("  Chi-square statistic:", round(chi_sq, 4), "\n")
cat("  Degrees of freedom: 2\n")
cat("  p-value:", round(p_value, 6), "\n")
cat("  Conclusion:", if(p_value > 0.05) 
  "Simulation matches analytical distribution" else 
    "Simulation differs significantly", "\n")

# Barplot comparison
barplot(comparison, beside = TRUE,
        col = c("skyblue", "salmon"),
        legend.text = rownames(comparison),
        args.legend = list(x = "topright"),
        main = "Traffic State Distribution at 6PM",
        xlab = "Traffic State", ylab = "Probability",
        ylim = c(0, max(comparison) + 0.05))
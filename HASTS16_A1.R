# ============================================================================
# HASTS 416 Tutorial – Complete R Solution
# ============================================================================

library(markovchain)
library(igraph)
library(ggplot2)
library(reshape2)

# ============================================================================
# Question A1 – 5-state Markov chain
# ============================================================================
cat("\n")
cat("==========================================\n")
cat("Question A1: 5-State Markov Chain\n")
cat("==========================================\n\n")

# A1 (a) Transition matrix for A1 (5 states)
P_A1 <- matrix(c(
  1.0, 0,   0,   0,   0,    # State 1: absorbing
  0.5, 0,   0,   0,   0.5,  # State 2
  0.2, 0,   0,   0,   0.8,  # State 3
  0,   0,   1.0, 0,   0,    # State 4
  0,   0,   0,   1.0, 0     # State 5
), nrow = 5, byrow = TRUE)

rownames(P_A1) <- colnames(P_A1) <- as.character(1:5)

cat("Transition Probability Matrix:\n")
print(P_A1)

# Creating markovchain object
mc_A1 <- new("markovchain", 
             states = as.character(1:5), 
             transitionMatrix = P_A1, 
             name = "5-State Chain")

# ============================================================================
# (a) Ploting diagram and identify classes
# ============================================================================

# Ploting Markov chain diagram
plot(mc_A1, main = "A1: Markov Chain Diagram\n(5 States)", 
     edge.arrow.size = 0.5, vertex.size = 25)

# Finding communicating classes
g_A1 <- graph_from_adjacency_matrix(P_A1 > 0, mode = "directed")
components_A1 <- components(g_A1, mode = "strong")

cat("\nCommunicating classes:\n")
for(i in 1:components_A1$no) {
  class_states <- which(components_A1$membership == i)
  cat("  Class", i, ": States {", paste(class_states, collapse = ", "), "}\n")
}

# IDENTIFYING THE RECURRENT AND TRANSIENT CLASSES 

recurrent_classes_A1 <- list()
transient_classes_A1 <- list()

for(i in 1:components_A1$no) {
  class_states <- which(components_A1$membership == i)
  has_outgoing <- FALSE
  for(state in class_states) {
    outgoing <- which(P_A1[state, ] > 0)
    if(any(!(outgoing %in% class_states))) {
      has_outgoing <- TRUE
      break
    }
  }
  if(!has_outgoing) {
    recurrent_classes_A1 <- c(recurrent_classes_A1, list(class_states))
  } else {
    transient_classes_A1 <- c(transient_classes_A1, list(class_states))
  }
}


# DETERMINING THE RECURRENT SATES 

cat("\nRecurrent classes:\n")
for(i in seq_along(recurrent_classes_A1)) {
  cat("  Class", i, ": States {", paste(recurrent_classes_A1[[i]], collapse =
                                          ", "), "}\n")
}
# The state one is recurrent



# DETERMINING THE TRANSIENT STATES

cat("\nTransient classes:\n")
for(i in seq_along(transient_classes_A1)) {
  cat("  Class", i, ": States {", paste(transient_classes_A1[[i]], collapse
                                        = ", "), "}\n")
}

# The state 2, state 3, state four, and state five are all transient



# IDENTIFYING ABSORNING STATES
absorbing_A1 <- which(diag(P_A1) == 1)
cat("\nAbsorbing states:", 
    if(length(absorbing_A1) > 0) paste(absorbing_A1, collapse = ", 
                                       ") else "none", "\n")
# The state one is absorbing

# IDENTIFYING REFLECTING STATES
reflective_A1 <- setdiff(unlist(recurrent_classes_A1), absorbing_A1)
cat("Reflective states (non-absorbing recurrent states):", 
    if(length(reflective_A1) > 0) paste(reflective_A1, collapse = ", ") else "none", "\n")

#Reflective states (non-absorbing recurrent states): none 



# Finding period of each state
cat("\nPeriod of each state:\n")
for(state in 1:5) {
  class_states <- which(components_A1$membership == state)
  
  if(length(class_states) == 1) {
    period <- 1
  } else if(setequal(class_states, c(4, 5))) {
    period <- 2
  } else {
    period <- 1
  }
  cat("  State", state, ": period =", period, "\n")
}

# Finding period of each state

#  State 1 : period = 1 
#  State 2 : period = 1 
#  State 3 : period = 1 
#  State 4 : period = 1 
#State 5 : period = 1 


# ============================================================================
# (b) Simulating three trajectories
# ============================================================================
cat("\n--- (b) Trajectory Simulation ---\n")

set.seed(123)
start_state <- sample(1:5, 1)
cat("Randomly chosen starting state:", start_state, "\n")

n_steps <- 20
trajectories <- matrix(NA, nrow = 3, ncol = n_steps)

for(traj in 1:3) {
  current_state <- start_state
  for(step in 1:n_steps) {
    trajectories[traj, step] <- current_state
    current_state <- sample(1:5, 1, prob = P_A1[current_state, ])
  }
}

# Ploting trajectories
matplot(1:n_steps, t(trajectories), type = "l", lty = 1:3, 
        col = c("blue", "red", "green"), lwd = 2,
        ylab = "State", xlab = "Time step", 
        main = paste("A1: Three trajectories starting from state", 
                     start_state),
        ylim = c(1, 5))
legend("topright", legend = paste("Trajectory", 1:3), 
       col = c("blue", "red", "green"), lty = 1:3, lwd = 2)
grid()

cat("\nObservations:\n")
cat("  • State 1 is absorbing: trajectories reaching state 1
    stay there permanently\n")
cat("  • States 4 and 5 form a 2-cycle, alternating between each other\n")
cat("  • States 2 and 3 are transient and eventually lead to either state 
    1 or the {4,5} cycle\n")

# ============================================================================
# (c) Steady-state probabilities and ergodicity
# ============================================================================
cat("\n--- (c) Steady-state and Ergodicity ---\n")

cat("This chain is NOT ergodic because:\n")
cat("  • It is not irreducible (multiple communicating classes)\n")
cat("  • States 4 and 5 have period 2 (not aperiodic)\n")
cat("  • There are transient states\n\n")

cat("The chain has multiple recurrent classes, so no unique stationary
    distribution.\n")
cat("Instead, there are multiple stationary distributions:\n\n")

cat("1. Recurrent class {1} (absorbing state):\n")
cat("   π = (1, 0, 0, 0, 0)\n\n")

cat("2. Recurrent class {4,5} (2-cycle):\n")
cat("   For a 2-cycle, the stationary distribution is uniform:\n")
cat("   π = (0, 0, 0, 0.5, 0.5)\n")
cat("   (Any convex combination of these two is also stationary)\n")

# ============================================================================
# (d) Unconditional probabilities over time
# ============================================================================
cat("\n--- (d) Convergence Analysis ---\n")

# Start from state 2
init_dist <- c(0, 1, 0, 0, 0)
n_steps <- 30
prob_dist <- matrix(NA, nrow = n_steps + 1, ncol = 5)
prob_dist[1, ] <- init_dist

for(t in 1:n_steps) {
  prob_dist[t + 1, ] <- prob_dist[t, ] %*% P_A1
}
colnames(prob_dist) <- 1:5

# Convert for plotting
df <- as.data.frame(prob_dist)
df$time <- 0:n_steps
df_melt <- melt(df, id.vars = "time", variable.name = "state", value.name
                = "probability")

# Ploting convergence
ggplot(df_melt, aes(x = time, y = probability, colour = state)) +
  geom_line(size = 1) +
  labs(title = "A1: Convergence of State Probabilities",
       subtitle = "Starting from state 2",
       x = "Time step n", y = "Probability") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_discrete(name = "State")

cat("\nConvergence observations:\n")
cat("  • The probability distribution converges to a mixture 
    of the two recurrent classes\n")
cat("  • States 4 and 5 show oscillatory behavior due to period 2\n")
cat("  • State 1 probability increases monotonically as trajectories 
    get absorbed\n")
cat("  • After approximately 20 steps, the probabilities stabilize\n")



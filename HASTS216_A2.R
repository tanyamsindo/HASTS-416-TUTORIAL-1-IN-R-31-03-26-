# ============================================================================
# HASTS 416 Tutorial – Complete R Solution
# ============================================================================

library(markovchain)
library(igraph)
library(ggplot2)
library(reshape2)


# ============================================================================
# Question A2 – 7-state Markov chain
# ============================================================================
cat("\n\n")
cat("==========================================\n")
cat("Question A2: 7-State Markov Chain\n")
cat("==========================================\n\n")

# Transition matrix for A2 (7 states)
P_A2 <- matrix(c(
  0,   1,   0,   0,   0,   0,   0,    # State 1
  1,   0,   0,   0,   0,   0,   0,    # State 2
  0,   0,   0,   0.4, 0.2, 0.2, 0.2,  # State 3
  0,   0,   0,   0,   0.2, 0.4, 0.4,  # State 4
  0.3, 0,   0,   0.1, 0.3, 0.1, 0.2,  # State 5
  0,   0,   0,   0.2, 0.2, 0.3, 0.3,  # State 6
  0,   0,   0,   0.5, 0.2, 0.2, 0.1   # State 7
), nrow = 7, byrow = TRUE)

rownames(P_A2) <- colnames(P_A2) <- as.character(1:7)

cat("Transition Probability Matrix:\n")
print(round(P_A2, 2))

# Create markovchain object
mc_A2 <- new("markovchain", 
             states = as.character(1:7), 
             transitionMatrix = P_A2, 
             name = "7-State Chain")

# ============================================================================
# (a) Plot diagram
# ============================================================================
cat("\n--- (a) Diagram ---\n")
plot(mc_A2, main = "A2: Markov Chain Diagram\n(7 States)", 
     edge.arrow.size = 0.4, vertex.size = 20)

# ============================================================================
# (b) Identifying classes, periods, absorbing/reflective states
# ============================================================================
cat("\n--- (b) Classification ---\n")

# Find communicating classes
g_A2 <- graph_from_adjacency_matrix(P_A2 > 0, mode = "directed")
components_A2 <- components(g_A2, mode = "strong")

cat("\nCommunicating classes:\n")
for(i in 1:components_A2$no) {
  class_states <- which(components_A2$membership == i)
  cat("  Class", i, ": States {", paste(class_states, collapse = ", "), "}\n")
}

#Class 1 : States { 3 }
#Class 2 : States { 4, 5, 6, 7 }
#Class 3 : States { 1, 2 }



# Identifying recurrent and transient classes
recurrent_classes_A2 <- list()
transient_classes_A2 <- list()

for(i in 1:components_A2$no) {
  class_states <- which(components_A2$membership == i)
  has_outgoing <- FALSE
  for(state in class_states) {
    outgoing <- which(P_A2[state, ] > 0)
    if(any(!(outgoing %in% class_states))) {
      has_outgoing <- TRUE
      break
    }
  }
  if(!has_outgoing) {
    recurrent_classes_A2 <- c(recurrent_classes_A2, list(class_states))
  } else {
    transient_classes_A2 <- c(transient_classes_A2, list(class_states))
  }
}

cat("\nRecurrent classes:\n")
for(i in seq_along(recurrent_classes_A2)) {
  cat("  Class", i, ": States {", paste(recurrent_classes_A2[[i]], collapse = ", "), "}\n")
}


# THE RECURRENT STATES
# Class 1 : States { 1, 2 }



cat("\nTransient classes:\n")
for(i in seq_along(transient_classes_A2)) {
  cat("  Class", i, ": States {", paste(transient_classes_A2[[i]], collapse = ", "), "}\n")
}

# THE TRANSIENT STATES


#Class 1 : States { 3 }
#Class 2 : States { 4, 5, 6, 7 }


# Identifying absorbing states
absorbing_A2 <- which(diag(P_A2) == 1)
cat("\nAbsorbing states:", 
    if(length(absorbing_A2) > 0) paste(absorbing_A2, collapse = ", ") 
    else "none", "\n")

# Absorbing states: none 


# Identify reflective states
reflective_A2 <- setdiff(unlist(recurrent_classes_A2), absorbing_A2)
cat("Reflective states (non-absorbing recurrent states):", 
    if(length(reflective_A2) > 0) paste(reflective_A2, collapse = ", ") 
    else "none", "\n")

# Reflective states (non-absorbing recurrent states): 1, 2 


# Finding period of each state
cat("\nPeriod of each state:\n")
for(state in 1:7) {
  class_states <- which(components_A2$membership == state)
  
  if(setequal(class_states, c(1, 2))) {
    period <- 2
  } else {
    period <- 1
  }
  cat("  State", state, ": period =", period, "\n")
}

cat("\nSummary:\n")
cat("  • States 1 and 2 form a recurrent class with period
    2 (2-cycle)\n")
cat("  • States 3-7 are all transient\n")
cat("  • No absorbing states (no state with self-loop
    probability 1)\n")
cat("  • No reflective states (the only recurrent
class has period 2,
    not reflective)\n")
cat("  • All transient states eventually lead to
    the {1,2} cycle\n")


#cat("\nSummary:\n")

#Summary:
# > cat("  • States 1 and 2 form a recurrent class with period
#    2 (2-cycle)\n")
#States 1 and 2 form a recurrent class with period
#2 (2-cycle)
#> cat("  • States 3-7 are all transient\n")
#• States 3-7 are all transient
#> cat("  • No absorbing states (no state with self-loop
#+     probability 1)\n")
#• No absorbing states (no state with self-loop
#                       probability 1)
#> cat("  • No reflective states (the only recurrent
#+ class has period 2,
#+     not reflective)\n")
#• No reflective states (the only recurrent
#                       class has period 2,
#                       not reflective)
#> cat("  • All transient states eventually lead to
#+     the {1,2} cycle\n")
#• All transient states eventually lead to
#the {1,2} cycle



# (c) Simulating two trajectories from a randomly selected state
# ============================================================================
cat("\n--- (c) Trajectory Simulation (A2) ---\n")

set.seed(456)
start_state_A2 <- sample(1:7, 1)
cat("Randomly chosen starting state:", start_state_A2, "\n")

n_steps_A2 <- 30
trajectories_A2 <- matrix(NA, nrow = 2, ncol = n_steps_A2)

for(traj in 1:2) {
  current_state <- start_state_A2
  for(step in 1:n_steps_A2) {
    trajectories_A2[traj, step] <- current_state
    current_state <- sample(1:7, 1, prob = P_A2[current_state, ])
  }
}

# Plot the two trajectories
matplot(1:n_steps_A2, t(trajectories_A2), type = "l", lty = 1:2,
        col = c("blue", "red"), lwd = 2,
        ylab = "State", xlab = "Time step",
        main = paste("A2: Two trajectories starting from state", start_state_A2),
        ylim = c(1, 7))
legend("topright", legend = paste("Trajectory", 1:2),
       col = c("blue", "red"), lty = 1:2, lwd = 2)
grid()

cat("\nObservations:\n")
cat("  • The chain eventually enters the recurrent class {1,2} (the 2-cycle).\n")
cat("  • Once in {1,2}, the chain alternates deterministically between states 1 and 2.\n")
cat("  • Transient states 3-7 may be visited initially, but after a few steps the chain\n")
cat("    gets trapped in the 2-cycle.\n")
cat("  • Because the recurrent class has period 2, the trajectories show a repeating pattern.\n")

# ============================================================================
# (d) Limiting probabilities and ergodicity
# ============================================================================
cat("\n--- (d) Limiting Probabilities and Ergodicity (A2) ---\n")

cat("This chain is NOT ergodic because:\n")
cat("  • It is not irreducible (multiple communicating classes)\n")
cat("  • The recurrent class {1,2} has period 2 (not aperiodic)\n")
cat("  • There are transient states\n\n")

cat("The limiting distribution does NOT exist in the usual sense because the chain is periodic.\n")
cat("However, we can compute the stationary distribution for the recurrent class:\n")

# Stationary distribution for the 2-cycle {1,2}
stationary_recurrent <- c(0.5, 0.5)
names(stationary_recurrent) <- c("1", "2")

cat("Stationary distribution for the recurrent class {1,2}:\n")
print(stationary_recurrent)

#> print(stationary_recurrent)
#1   2 
#0.5 0.5 


cat("\nInterpretation:\n")
cat("  • If we only look at the times when the chain is in the recurrent class,\n")
cat("    the long-run proportion of time spent in state 1 is 0.5, and in state 2 is 0.5.\n")
cat("  • Because the chain is periodic, the probabilities do not converge; instead they oscillate.\n")
cat("  • Starting from a state in the recurrent class, the distribution at time n is:\n")
cat("      - If n is even: probability 1 on state 1 (if started at state 1) or state 2 (if started at state 2)\n")
cat("      - If n is odd: probability 1 on the opposite state.\n")
cat("  • Starting from a transient state, the probabilities eventually become concentrated on\n")
cat("    the recurrent class, but still oscillate between states 1 and 2.\n")

# Absorption probabilities
cat("\nAbsorption probabilities into the recurrent class {1,2} from each transient state:\n")
transient_idx <- c(3,4,5,6,7)
recurrent_idx <- c(1,2)
Q <- P_A2[transient_idx, transient_idx]
R <- P_A2[transient_idx, recurrent_idx]
I <- diag(length(transient_idx))
N <- solve(I - Q)
B <- N %*% R
rownames(B) <- transient_idx
colnames(B) <- recurrent_idx
print(round(B, 4))

cat("\nExpected time until absorption into {1,2}:\n")
expected_time <- rowSums(N)
names(expected_time) <- transient_idx
print(round(expected_time, 2))


#> print(round(expected_time, 2))
#3  4  5  6  7 
#15 15 10 15 15 





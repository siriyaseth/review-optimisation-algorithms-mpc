# Review of Optimisation Algorithms for MPC
This repository contains the code for my final project in Introduction to Optimisation.
Three algorithms are implemented in MATLAB and applied to some examples:
- Projected Gradient
- Dual Ascent with Decomposition
- Augmented Lagrangian

## Running Examples
The following files generate the figures inside the final report.

For gradient projection applied to MPC with the ball-and-plate system:
- gradient-projection/simulate_closed-loop.m - Closed-loop simulation with state and input figures
- gradient-projection/test_warmstart.m - Number of iterations before convergence threshold satisfied over closed-loop simulation with and without warm starting
- gradient-projection/test_convergence.m - Number of iterations before convergence threshold satisfied for varying step sizes

For dual decomposition applied to MPC with the distributed chicken coop control:
- dual-decomposition/simulate_closed-loop.m - Closed-loop simulation with state and input figures
- dual-decomposition/test-chicken_coop-warmstart - Number of iterations before convergence threshold satisfied over closed-loop simulation with and without warm starting

For augmented lagrangian applied to MPC with the CSTR control:
- augmented-lagrangian/test_cstr_auglag.m - Closed-loop simulation with state and input figures as well as iterations before convergence threshold satisfied

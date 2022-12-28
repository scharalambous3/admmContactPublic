# Alternating Direction Method of Multipliers (ADMM)
Extends prior work, *Aydinoglu, Alp, and Michael Posa. "Real-time multi-contact model predictive control via admm." 2022 International Conference on Robotics and Automation (ICRA). IEEE, 2022.*, for use of ADMM for multi-contact systems. Currently implemented:
(1) finger gaiting, (2) cart-pole with soft walls, (3) 2D floating base, and (4) an alternative formulation of 2D floating base

## Features
- Problem agnostic. Can extend to any new system by defining a parameters file
- Implicit time-stepping is used for modelling (1) and (4) 
- Calculation of primal and dual residuals for monitoring convergence

## Results
*These results are obtained using modifications to the original work. These modifications are currently private*
<p style="text-align:center;"><img src="https://github.com/scharalambous3/admmContactPublic/blob/main/results.png" alt="Logo"></p>

**Functions**:
- projectionSubproblem: Project solution of QP onto the linear complementarity constraints
- projectionSubproblemGrouping Project solution of QP onto the linear complementarity constraints. Includes grouping of integer decision variables
- benchmarks/: Benchmarks to compare ADMM results against, using MINLP or NLP formulations
- solveADMM: Iterate through convex subproblem, projection and dual update until convergence or max iterations
- solveConvexSubproblem: Solve convex subproblem (QP since linear complementarity constraints are not included here)
- plotPerf: Plots metrics (orthogonality violation, objective value of given state trajectory, primal and dual residuals), trajectory of state variables and boolean assignment of the integer decision variables
- plotSeparation: Plots variables involved in separation complementarity constraint
- plotStates: Plots state variables involved in complementarity constraints
- params/: Parameters for each system of the different problem types
- util/: Utility functions for rollout, visualization, discrete dynamics transition map, 2D cross-product and calculation of violation



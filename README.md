# zika-vaccination
Fortran 90 code for the work entitled "Optimal control of vaccination in a vector-borne reaction-diffusion model for Zika virus".

The code solves the optimal control PDE problem using:
- Finite element method;
- Finite difference (implicit Euler) method;
- Predictor-corrector method for successive linearization in time stepping;
- Forward-backward sweep for the iterative optimal solution;
- Conjugate gradient with diagonal preconditioner for linear system solutions with sparse storage;
- BOBYQA package (Powell, 2009) for optimization in the non linear regression of parameters.

Mesh files are necessary for a particular finite element mesh.

The code can also provide Latin Hypercube Sampling solutions used to perform Partial Rank Correlation Coefficient analysis.

# zika-vaccination
Fortran code for the work entitled "Optimal control of vaccination in a vector-borne reaction-diffusion model for Zika virus".

The code solves the optimal control PDE problem using:
- Finite element method;
- Finite difference (implicit Euler) method;
- Predictor-corrector method for successive linearization in time stepping;
- Forward-backward sweep for the iterative optimal solution;
- BOBYQA package (Powell, 2009) for optimization in the non linear regression of parameters.

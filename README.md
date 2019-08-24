# zika-vaccination
Fortran 90 code for the work entitled "Optimal control of vaccination in a vector-borne reaction-diffusion model for Zika virus".

This code does not receive maintenance and is available here only for reference, use at own risk.

The code solves the optimal control PDE problem using:
- Finite element method;
- Finite difference (implicit Euler) method;
- Predictor-corrector method for successive linearization in time stepping;
- Forward-backward sweep for the iterative optimal solution;
- Conjugate gradient with diagonal preconditioner for linear system solutions with sparse storage;
- BOBYQA package (Powell, 2009) for optimization in the non linear regression of parameters.

Mesh files are necessary for a particular finite element mesh.

The code can also provide Latin Hypercube Sampling solutions used to perform Partial Rank Correlation Coefficient analysis.

References:
- http://repositorio.unicamp.br/jspui/handle/REPOSIP/334344?mode=full (PhD. thesis)
- https://link.springer.com/article/10.1007/s00285-019-01390-z

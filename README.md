# GPspline
(Experimental) Package for a varying-coefficients spline regression with GP-regularization of the spline coefficients using a separate set of covariates

The intended use case is a (for now) single dimensional set of continuous variable(s) Z whose marginal (predictive/causal) effect we are interested. The other set, the control variables X or "confounders" for causal inference, are used to make the function approximation more accurate and for the causal effects control for confounding. Using Gaussian processes, we can use differentiable spline bases to obtain the marginal effect.

Formally,
```
y = g(x) * b(z) + eps
```
where ```g``` has a GP-prior and ```b``` is a polynomial spline design vector.

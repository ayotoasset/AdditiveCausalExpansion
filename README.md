# GPspline
(Experimental) Package for a varying-coefficients spline regression with GP-regularization of the spline coefficients using a separate set of covariates. The (GP) hyperparameters are chosen using empirical Bayes and note that the spline knots are fixed in number and location.

The intended use case is a (for now) single dimensional set of continuous variable(s) Z whose marginal (predictive/causal) effect we are interested. The other set, the control variables X or "confounders" for causal inference, are used to make the function approximation more accurate and for the causal effects control for confounding. Using Gaussian processes, we can use differentiable spline bases to obtain the marginal effect.

Formally,
```
y = mu + g(x) * b(z) + eps,
```
where each element of ```g``` has an independent GP-prior with covariance kernel ```K_g``` and zero prior-mean, ```b``` is a polynomial spline design vector (with dimension ```B```), and ```mu``` is the mean. The noise term ```eps``` is Gaussian with unknown variance ```sig^2```.

We can write the model in reduced form ```y = f(x,z) + eps``` with ```f ~ GP(mu, K_r)``` where the additive kernel is given by
```
K_r(i,j) = sum_{l=1}^B K_{g_l}(x_i,x_j) b_l(z_i) b_l(z_j).
```
This constitutes a proper kernel (sum of a product of kernels) and we can use standard Gaussian process inference methods to obtain the posterior distribution using empirical Bayes.

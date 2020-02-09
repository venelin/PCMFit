# Copied from this thread on stats.stackexchange:

# https://stats.stackexchange.com/questions/107597/is-there-a-way-to-use-the-covariance-matrix-to-find-coefficients-for-multiple-re

#
# 1. Generate some data.
#
n <- 10        # Data set size
p <- 3         # Number of regressors
set.seed(17)
z <- matrix(rnorm(n*(p+1)), nrow=n, dimnames=list(NULL, paste0("x", 1:(p+1))))
y <- z[, p+1]
x <- z[, -(p+1), drop=FALSE];
#
# 2. Find the OLS coefficients from the covariances only.
#
a <- cov(x)
b <- cov(x,y)
beta.hat <- solve(a, b)[, 1]  # Coefficients from the covariance matrix
#
# 2a. Find the intercept from the means and coefficients.
#
y.bar <- mean(y)
x.bar <- colMeans(x)
intercept <- y.bar - x.bar %*% beta.hat


(rbind(`From covariances` = c(`(Intercept)`=intercept, beta.hat),
       `From data via OLS` = coef(lm(y ~ x))))

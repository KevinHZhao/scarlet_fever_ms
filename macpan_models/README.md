fit.R: fits logistic rate parameter, uses sparse matrix multiplication, but no random effects (due to memory usage), penalizes the fourier coefficient (for sin(nx), cos(nx), and beta_0)

fit_nopenalize.R: same as fit.R but no penalization, much fewer fit parameters but displays overfitting

fit_sensitivity.R: fit.R but tests with different conditions for CFP parameters (81 different combinations, see manuscript initial conditions section in appendix), now using iterations to see what a get starting value is.

Right now, we use a few iterations of fit_sensitivity.R to get the best optimum that we can find (-loglik = ~29k) and then run `finalize_sensitivity.R` to move the best fit to the `ms_data` folder.

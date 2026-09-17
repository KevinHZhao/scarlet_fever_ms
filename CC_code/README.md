fit.R: fits logistic rate parameter, uses sparse matrix multiplication, but no random effects (due to memory usage), penalizes the fourier coefficient (for sin(nx), cos(nx), and beta_0)

fit_nopenalize.R: same as fit.R but no penalization, much fewer fit parameters but displays overfitting

fit_sensitivity.R: fit.R but tests with different conditions for CFP parameters (81 different combinations, see manuscript initial conditions section in appendix)

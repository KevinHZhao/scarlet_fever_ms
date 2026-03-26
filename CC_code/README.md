fit.R: fits logistic rate parameter, uses sparse matrix multiplication, but no random effects (due to memory usage), penalizes the fourier coefficient (for sin(nx), cos(nx), and beta_0)

fit_nopenalize.R: same as fit.R but no penalization, much fewer fit parameters but displays overfitting

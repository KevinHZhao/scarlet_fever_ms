# Transition analysis for scarlet fever paper

1.  First, run `fit.R` and `fit_nopenalize.R` in `CC_code/` to generate the `macpan2` model output using RGWR scarlet fever data (output will be in `CC_code/output`).
2.  Then, we need to determine the breaks for `breaks.csv`, we do this using the `boxplots.R` script in `breaks/`, and requires the user to make a decision qualitatively for which break points make the most sense.
3.  Once `breaks.csv` is ready, use XPPAUT to use the method of Papst and Earn to do alpha continuation, see `alpha_continuation/steps.md`.
4.  Afterwards, run the relevant ODE simulations for attractor and transient periods in `TwoParSim/`, see `TwoParSim/README.md`.
5.  To run the `fastbeta` comparison between the `macpan2` model and the equivalent `fastbeta` model, run the scripts in `fastbeta_compare/`, see `fastbeta_compare/README.md`.

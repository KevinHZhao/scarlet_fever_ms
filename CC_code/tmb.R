library(macpan2)

computations = list(
  N ~ pop[time_step(1)],
  model_pop ~ sum(S, I, R),
  outflows ~ model_pop - N ## term accounting for deaths and immigration/emigration (can think of as a reduction to death rate, might be positive)
)

flows = list(
  mp_per_capita_flow("S", "I", "I * beta / model_pop", "infection")
  , mp_per_capita_flow("I", "R", recovery ~ (1 - CFP) * gamma)
  , mp_per_capita_flow("I", "D", death ~ CFP * gamma)
  , mp_per_capita_outflow("S", mu_S ~ (outflows)/model_pop)
  , mp_per_capita_outflow("I", mu_I ~ (outflows)/model_pop)
  , mp_per_capita_outflow("R", mu_R ~ (outflows)/model_pop)
  , mp_inflow("S", nu ~ births[time_step(1)])
) # Use mu_S and mu_I and mu_R for outflows, can't directly obtain (can look at mp_expand()), can trust that nu is exactly births

## set defaults
default = list(
  beta = 0.2
  , gamma = 7/15
  , births = 0
  #, acm = 0
  , pop = 0
  , outflows = 0
  , N = 100
  , S = 99
  , I = 1
  , R = 0
  , P = 0
  , D = 0
  , CFP = 0.01
)

## model specification
spec = mp_tmb_model_spec(
  during = c(computations,flows)
  , default = default
)

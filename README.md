# random_effect_ODE

main.R: code for analyzing the duck flu data based on Bayesian multilevel mixed-effects model. The outputs include the posterior distribution for the strain level random effect. The fitted curve for the virus load of each individual is also included.  

ode_mcmc.R: code for updating ODE variable and parameter for one individual

ode_fit.R: code for finding the optimal ODE fit for the virus load of each individual 

ode_formula.R: code for providing ODE formula and its contribution to the update of spline coefficients in MCMC procedure 

camille-duck-data.csv: Duck flu data file which includes 6 strains and each strain includes 4 or 5 individuals

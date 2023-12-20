## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi = 72) # set global knitr options

## -----------------------------------------------------------------------------
library(remstimate)

## -----------------------------------------------------------------------------
# setting `ncores` to 1 (the user can change this parameter)
ncores <- 1L

# loading data
data(tie_data)

# true parameters' values
tie_data$true.pars

# processing the event sequence with 'remify'
tie_reh <- remify::remify(edgelist = tie_data$edgelist, model = "tie")

# summary of the (processed) relational event network
summary(tie_reh)

## -----------------------------------------------------------------------------
# specifying linear predictor (with `remstats`) using a 'formula'
tie_model <- ~ 1 + remstats::indegreeSender() + 
              remstats::inertia() + remstats::reciprocity() 

## -----------------------------------------------------------------------------
# calculating statistics (with `remstats`)
tie_stats <- remstats::remstats(reh = tie_reh, tie_effects = tie_model)

# the 'tomstats' 'remstats' object
tie_stats


## -----------------------------------------------------------------------------
# for example the method "MLE"
remstimate::remstimate(reh = tie_reh,
                          stats =  tie_stats,
                          method = "MLE",
                          ncores = ncores)    

## -----------------------------------------------------------------------------
  tie_mle <- remstimate::remstimate(reh = tie_reh,
                          stats = tie_stats,
                          ncores = ncores,
                          method = "MLE",
                          WAIC = TRUE, # setting WAIC computation to TRUE
                          nsimWAIC = 100) # number of draws for the computation of the WAIC set to 100         

## -----------------------------------------------------------------------------
# printing the 'remstimate' object
tie_mle

## -----------------------------------------------------------------------------
# summary of the 'remstimate' object
summary(tie_mle)

## -----------------------------------------------------------------------------
# aic
aic(tie_mle)

# aicc
aicc(tie_mle)
  
# bic 
bic(tie_mle)

#waic 
waic(tie_mle)

## -----------------------------------------------------------------------------
# diagnostics
tie_mle_diagnostics <- diagnostics(object = tie_mle, reh = tie_reh, stats = tie_stats)

## ----out.width="50%", dev=c("jpeg")-------------------------------------------
# plot
plot(x = tie_mle, reh  = tie_reh, diagnostics = tie_mle_diagnostics)

## -----------------------------------------------------------------------------
tie_gd <- remstimate::remstimate(reh = tie_reh,
                        stats =  tie_stats,
                        ncores = ncores,
                        method = "GDADAMAX",
                        epochs = 200L, # number of iterations for the Gradient-Descent algorithm
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100) # number of draws for the computation of the WAIC set to 100     
# print 
tie_gd

## -----------------------------------------------------------------------------
# diagnostics
tie_gd_diagnostics <- diagnostics(object = tie_gd, reh = tie_reh, stats = tie_stats)
# plot
# plot(x = tie_gd, reh  = tie_reh, diagnostics = tie_gd_diagnostics) # uncomment to use the plot function

## -----------------------------------------------------------------------------
library(mvnfast) # loading package for fast simulation from a multivariate Student t distribution
priormvt <- mvnfast::dmvt # defining which distribution we want to use from the 'mvnfast' package
tie_bsir <- remstimate::remstimate(reh = tie_reh,
                        stats =  tie_stats,
                        ncores = ncores,
                        method = "BSIR",
                        nsim = 200L, # 200 draws from the posterior distribution
                        prior = priormvt, # defining prior here, prior parameters follow below
                        mu = rep(0,dim(tie_stats)[3]), # prior mu value
                        sigma = diag(dim(tie_stats)[3])*1.5, # prior sigma value
                        df = 1, # prior df value
                        log = TRUE, # requiring log density values from the prior,
                        seed = 23029, # set a seed only for reproducibility purposes
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100 # number of draws for the computation of the WAIC set to 100     
                        )

# summary 
summary(tie_bsir)

## -----------------------------------------------------------------------------
# diagnostics
tie_bsir_diagnostics <- diagnostics(object = tie_bsir, reh = tie_reh, stats = tie_stats)
# plot
# plot(x = tie_bsir, reh  = tie_reh, diagnostics = tie_bsir_diagnostics) # uncomment to use the plot function

## -----------------------------------------------------------------------------
tie_hmc <- remstimate::remstimate(reh = tie_reh,
                        stats =  tie_stats,
                        method = "HMC",
                        ncores = ncores,
                        nsim = 200L, # 200 draws to generate per each chain
                        nchains = 4L, # 4 chains to generate
                        burnin = 200L, # burnin length is 200
                        thin = 2L, # thinning size set to 2 (the final length of the chains will be 100)
                        seed = 23029, # set a seed only for reproducibility purposes
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100 # number of draws for the computation of the WAIC set to 100     
                        )

# summary 
summary(tie_hmc)

## ----out.width="50%", dev=c("jpeg")-------------------------------------------
# diagnostics
tie_hmc_diagnostics <- diagnostics(object = tie_hmc, reh = tie_reh, stats = tie_stats)
# plot (histograms and trace plot have highest posterior density intervals dashed lines in blue and posterior estimate in red)
plot(x = tie_hmc, reh  = tie_reh, diagnostics = tie_hmc_diagnostics)

## -----------------------------------------------------------------------------
# setting `ncores` to 1 (the user can change this parameter)
ncores <- 1L

# loading data
data(ao_data)

# true parameters' values
ao_data$true.pars

# processing event sequence with 'remify'
ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")

# summary of the relational event network
summary(ao_reh)

## -----------------------------------------------------------------------------
# specifying linear predictor (for rate and choice model, with `remstats`)
rate_model <- ~ 1 + remstats::indegreeSender()
choice_model <- ~ remstats::inertia() + remstats::reciprocity()


## -----------------------------------------------------------------------------
# calculating statistics (with `remstats`)
ao_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model)

# the 'aomstats' 'remstats' object
ao_stats

## -----------------------------------------------------------------------------
# for example the method "MLE"
remstimate::remstimate(reh = ao_reh,
                          stats =  ao_stats,
                          method = "MLE",
                          ncores = ncores)    

## -----------------------------------------------------------------------------
ao_mle <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_stats,
                        ncores = ncores,
                        method = "MLE",
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100) # number of draws for the computation of the WAIC set to 100            

## -----------------------------------------------------------------------------
# printing the 'remstimate' object 
ao_mle

## -----------------------------------------------------------------------------
# summary of the 'remstimate' object
summary(ao_mle)

## -----------------------------------------------------------------------------
# aic
aic(ao_mle)

# aicc
aicc(ao_mle)
  
# bic 
bic(ao_mle)

#waic 
waic(ao_mle)

## -----------------------------------------------------------------------------
# diagnostics
ao_mle_diagnostics <- diagnostics(object = ao_mle, reh = ao_reh, stats = ao_stats)

## ----out.width="50%", dev=c("jpeg")-------------------------------------------
# plot
plot(x = ao_mle, reh  = ao_reh, diagnostics = ao_mle_diagnostics)

## -----------------------------------------------------------------------------
ao_gd <- remstimate::remstimate(reh = ao_reh,
                        stats =  ao_stats,
                        ncores = ncores,
                        method = "GDADAMAX",
                        epochs = 200L, # number of iterations of the Gradient-Descent algorithm
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100) # number of draws for the computation of the WAIC set to 100     
# print 
ao_gd

## -----------------------------------------------------------------------------
# diagnostics
ao_gd_diagnostics <- diagnostics(object = ao_gd, reh = ao_reh, stats = ao_stats)
# plot
# plot(x = ao_gd, reh  = ao_reh, diagnostics = ao_gd_diagnostics) # uncomment to use the plot function

## -----------------------------------------------------------------------------
library(mvnfast) # loading package for fast simulation from a multivariate Student t distribution
priormvt <- mvnfast::dmvt # defining which distribution we want to use from the 'mvnfast' package
ao_bsir <- remstimate::remstimate(reh = ao_reh,
                        stats =  ao_stats,
                        ncores = ncores,
                        method = "BSIR",
                        nsim = 100L, # 100 draws from the posterior distribution
                        prior = list(sender_model = priormvt, receiver_model = priormvt), #  defining prior here, prior parameters follow below
                        prior_args = list(sender_model =  list(mu = rep(0,dim(ao_stats$sender_stats)[3]), # prior mu value for sender_model
                                                            sigma = diag(dim(ao_stats$sender_stats)[3])*1.5, # prior sigma value for sender_model
                                                            df = 1),  # prior df value
                                        receiver_model = list(mu = rep(0,dim(ao_stats$receiver_stats)[3]), # prior mu value for receiver_model
                                                            sigma = diag(dim(ao_stats$receiver_stats)[3])*1.5, # prior sigma value for receiver_model
                                                            df = 1)), # prior df value
                        log = TRUE, # requiring log density values from the prior,
                        seed = 20929, # set a seed only for reproducibility purposes
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100 # number of draws for the computation of the WAIC set to 100     
                        )

# summary 
summary(ao_bsir)

## -----------------------------------------------------------------------------
# diagnostics
ao_bsir_diagnostics <- diagnostics(object = ao_bsir, reh = ao_reh, stats = ao_stats)
# plot (only for the receiver_model, by setting sender_model = NA)
# plot(x = ao_bsir, reh  = ao_reh, diagnostics = ao_bsir_diagnostics, sender_model = NA) # uncomment to use the plot function

## -----------------------------------------------------------------------------
ao_hmc <- remstimate::remstimate(reh = ao_reh,
                        stats =  ao_stats,
                        method = "HMC",
                        ncores = ncores,
                        nsim = 300L, # 300 draws to generate per each chain
                        nchains = 4L, # 4 chains (each one long 200 draws) to generate
                        burnin = 300L, # burnin length is 300
                        L = 100L, # number of leap-frog steps
                        epsilon = 0.1/100, # size of a leap-frog step
                        thin = 2L, # thinning size (this will reduce the final length of each chain will be 150)
                        seed = 23029, # set a seed only for reproducibility purposes
                        WAIC = TRUE, # setting WAIC computation to TRUE
                        nsimWAIC = 100 # number of draws for the computation of the WAIC set to 100     
                        )

# summary 
summary(ao_hmc)

## ----out.width="50%", dev=c("jpeg")-------------------------------------------
# diagnostics
ao_hmc_diagnostics <- diagnostics(object = ao_hmc, reh = ao_reh, stats = ao_stats)
# plot (only for the receiver_model, by setting sender_model = NA)
plot(x = ao_hmc, reh  = ao_reh, diagnostics = ao_hmc_diagnostics, sender_model = NA)


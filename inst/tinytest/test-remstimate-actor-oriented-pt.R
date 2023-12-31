## testing actor-oriented modeling ##

# loading data
data(ao_data)

set_seed <- 23929

# processing data with simultaneous events (two events observed at the same time point)
ao_data$edgelist$time <- floor(ao_data$edgelist$time) 
ao_data$edgelist$time[seq(5,95,by=5)] <- ao_data$edgelist$time[seq(5,95,by=5)-1]
ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor")

# specifying linear predictor (for rate and choice model)
rate_model <- ~ 1 + remstats::indegreeSender()
choice_model <- ~ remstats::inertia() + remstats::reciprocity()

# calculating statistics
ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model, method="pt")

# (1) method = "MLE"
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "MLE"))                      
ao_mle <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "MLE")
expect_silent(ao_mle)
expect_inherits(ao_mle,"remstimate")  
expect_length(ao_mle,2)
expect_identical(names(ao_mle),c("sender_model","receiver_model"))
expect_length(ao_mle$sender_model,17)
expect_length(ao_mle$receiver_model,17)
expect_identical(names(ao_mle$sender_model),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_identical(names(ao_mle$receiver_model),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_length(attributes(ao_mle),10)
expect_identical(names(attributes(ao_mle)),c("names","class","formula","model","ordinal","method","approach","statistics","where_is_baseline","ncores"))
expect_identical(attr(ao_mle,"approach"),"Frequentist")
expect_silent(print(ao_mle))
expect_silent(summary(ao_mle))
expect_silent(diagnostics(object = ao_mle, reh = ao_reh, stats = ao_reh_stats))
ao_reh_diagnostics <- diagnostics(object = ao_mle, reh = ao_reh, stats = ao_reh_stats)
expect_silent(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics))
expect_silent(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics, which = 2, sender_effects = "indegreeSender", receiver_effects = "inertia")) # plotting specific effects from both models
expect_silent(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics, which = 2, sender_effects = NA, receiver_effects = "inertia")) # plotting specific effect from one model only 

expect_silent(aic(ao_mle))
expect_silent(aicc(ao_mle))
expect_silent(bic(ao_mle))
expect_silent(waic(ao_mle))

# error when 'diagnostics' has at least one statistic's name wrong
colnames(ao_reh_diagnostics$sender_model$residuals$smoothing_weights)[1] <- "INDEGREEsender"
expect_error(plot(x = ao_mle,reh = ao_reh, diagnostics = ao_reh_diagnostics),
"one or more effects not found inside the object 'diagnostics'.",
fixed=TRUE)

# test diagnostics with only baseline
ao_stats_baseline <- remstats::remstats(reh = ao_reh, sender_effects = ~1 , method = "pt")
ao_mle_only_baseline <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_stats_baseline,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = ao_mle_only_baseline, reh = ao_reh, stats = ao_stats_baseline))



# tests on "WAIC" for "MLE"

# WAIC = TRUE
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE))
# nsimWAIC is not a numeric scalar                       
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = "text"))
# nsimWAIC is supplied                        
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
ao_mle_with_waic <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100)
expect_silent(print(ao_mle_with_waic))
expect_silent(summary(ao_mle_with_waic))    
expect_silent(waic(ao_mle_with_waic))      

# WAIC for ordinal likelihod + testing omit_dyad routine for ordinal likelihood and tie-oriented model
omit_dyad <- list()
omit_dyad[[1]] <- list(time = c(330,585), dyad = data.frame(actor1=c(NA,"4"),actor2=c("4",NA),type=c(NA,NA))) # excluding actor 4 from the risk set between (observed) time points 330 and 585
ao_reh_ordinal <- remify::remify(edgelist = ao_data$edgelist, model = "actor", ordinal = TRUE, omit_dyad = omit_dyad)
ao_reh_stats_ordinal <- remstats::remstats(reh = ao_reh_ordinal, sender_effects = rate_model, receiver_effects = choice_model, method="pt")
expect_silent(remstimate::remstimate(reh = ao_reh_ordinal,
                        stats = ao_reh_stats_ordinal,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
#testing only ordinal likelihood for actor-oriented model
ao_mle_ordinal_omit_dyad <- remstimate::remstimate(reh = ao_reh_ordinal,
                        stats = ao_reh_stats_ordinal,
                        ncores = 1L)                     
# WAIC for interval likelihood with omit_dyad
ao_reh_omit <- remify::remify(edgelist = ao_data$edgelist, model = "actor", omit_dyad = omit_dyad)           
ao_reh_stats_omit <- remstats::remstats(reh = ao_reh_omit, sender_effects = rate_model, receiver_effects = choice_model, method="pt")
expect_silent(remstimate::remstimate(reh = ao_reh_omit,
                        stats = ao_reh_stats_omit,
                        ncores = 1L,
                        WAIC = TRUE,
                        nsimWAIC = 100))
expect_silent(diagnostics(object = ao_mle_ordinal_omit_dyad,reh = ao_reh_ordinal,stats = ao_reh_stats_ordinal))
# (2) method = "GDADAMAX" 
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10,
                        epsilon = NULL))
ao_gdadamax <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10,
                        epsilon = NULL)
expect_silent(ao_gdadamax)
expect_inherits(ao_gdadamax,"remstimate")  
expect_length(ao_gdadamax,2)
expect_identical(names(ao_gdadamax),c("sender_model","receiver_model"))
# (2) method = "GDADAMAX" 
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10,
                        epsilon = NULL))
ao_gdadamax <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10,
                        epsilon = NULL)
expect_silent(ao_gdadamax)
expect_inherits(ao_gdadamax,"remstimate")  
expect_length(ao_gdadamax,2)
expect_identical(names(ao_gdadamax),c("sender_model","receiver_model"))
expect_length(ao_gdadamax$sender_model,17)
expect_length(ao_gdadamax$receiver_model,17)
expect_identical(names(ao_gdadamax$sender_model),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_identical(names(ao_gdadamax$receiver_model),c("coefficients","loglik","gradient","hessian","vcov","se","residual.deviance","null.deviance","model.deviance","df.null","df.model","df.residual","AIC","AICC","BIC","converged","iterations"))
expect_length(attributes(ao_gdadamax),13)
expect_identical(names(attributes(ao_gdadamax)),c("names","class","formula","model","ordinal","method","approach","statistics","epochs","epsilon","init","where_is_baseline","ncores"))
expect_identical(attr(ao_gdadamax,"approach"),"Frequentist")
expect_silent(print(ao_gdadamax))
expect_silent(summary(ao_gdadamax))
expect_silent(diagnostics(object = ao_gdadamax, reh = ao_reh, stats = ao_reh_stats))
ao_reh_diagnostics <- diagnostics(object = ao_gdadamax, reh = ao_reh, stats = ao_reh_stats)
expect_silent(plot(x = ao_gdadamax,reh = ao_reh, diagnostics = ao_reh_diagnostics))
expect_silent(aic(ao_gdadamax))
expect_silent(aicc(ao_gdadamax))
expect_silent(bic(ao_gdadamax))
expect_silent(waic(ao_gdadamax)) 

# (3) method  = "BSIR" 

## (3.1) with prior = NULL
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL,
                        seed = set_seed))
ao_bsir_no_prior <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL,
                        seed = set_seed)
expect_silent(ao_bsir_no_prior)
expect_inherits(ao_bsir_no_prior,"remstimate")  
expect_length(ao_bsir_no_prior,2)
expect_identical(names(ao_bsir_no_prior),c("sender_model","receiver_model"))
expect_length(ao_bsir_no_prior$sender_model,12)
expect_length(ao_bsir_no_prior$receiver_model,12)
expect_identical(names(ao_bsir_no_prior$sender_model),c("log_posterior","draws","log_proposal","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual"))
expect_identical(names(ao_bsir_no_prior$receiver_model),c("log_posterior","draws","log_proposal","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual"))
expect_length(attributes(ao_bsir_no_prior),12)
expect_identical(names(attributes(ao_bsir_no_prior)),c("names","class","formula","model","ordinal","method","approach","statistics","nsim","seed","where_is_baseline","ncores"))
expect_identical(attr(ao_bsir_no_prior,"approach"),"Bayesian")
expect_silent(print(ao_bsir_no_prior))
expect_silent(summary(ao_bsir_no_prior))
expect_silent(diagnostics(object = ao_bsir_no_prior, reh = ao_reh, stats = ao_reh_stats))
ao_reh_diagnostics <- diagnostics(object = ao_bsir_no_prior, reh = ao_reh, stats = ao_reh_stats)
expect_silent(plot(x = ao_bsir_no_prior,reh = ao_reh, diagnostics = ao_reh_diagnostics))
expect_error(aic(ao_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(ao_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(ao_bsir_no_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_silent(waic(ao_bsir_no_prior)) 

# tests on "WAIC" for "BSIR"                       
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 100L,
                        prior = NULL,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10))
ao_reh_bsir_waic <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 100L,
                        prior = NULL,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10)
expect_silent(summary(ao_reh_bsir_waic))

## (3.2) with a specified prior
prior_list <- list(sender_model = mvnfast::dmvt, receiver_model = mvnfast::dmvn)
prior_args <- list(sender_model = list(mu = rep(0,dim(ao_reh_stats$sender_stats)[3]), sigma = diag(dim(ao_reh_stats$sender_stats)[3]), df = 1),
 receiver_model = list(mu = rep(0,dim(ao_reh_stats$sender_stats)[3]), sigma = diag(dim(ao_reh_stats$receiver_stats)[3])))
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = prior_list,
                        prior_args = prior_args,
                        seed = set_seed))
ao_bsir_with_prior <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = prior_list,
                        prior_args = prior_args,
                        seed = set_seed)
expect_silent(ao_bsir_with_prior)
expect_inherits(ao_bsir_with_prior,"remstimate")  
expect_length(ao_bsir_with_prior,2)
expect_identical(names(ao_bsir_with_prior),c("sender_model","receiver_model"))
expect_length(ao_bsir_with_prior$sender_model,13)
expect_length(ao_bsir_with_prior$receiver_model,13)
expect_identical(names(ao_bsir_with_prior$sender_model),c("log_posterior","draws","log_proposal","log_prior","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual"))
expect_identical(names(ao_bsir_with_prior$receiver_model),c("log_posterior","draws","log_proposal","log_prior","irw","coefficients","loglik","post.mean","vcov","sd","df.null","df.model","df.residual"))
expect_length(attributes(ao_bsir_with_prior),13)
expect_identical(names(attributes(ao_bsir_with_prior)),c("names","class","formula","model","ordinal","method","approach","statistics","prior","nsim","seed","where_is_baseline","ncores"))
expect_identical(attr(ao_bsir_with_prior,"approach"),"Bayesian")
expect_silent(print(ao_bsir_with_prior))
expect_silent(summary(ao_bsir_with_prior))
expect_silent(diagnostics(object = ao_bsir_with_prior, reh = ao_reh, stats = ao_reh_stats))
ao_reh_diagnostics <- diagnostics(object = ao_bsir_with_prior, reh = ao_reh, stats = ao_reh_stats)
expect_silent(plot(x = ao_bsir_with_prior,reh = ao_reh, diagnostics = ao_reh_diagnostics))
expect_error(aic(ao_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(ao_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(ao_bsir_with_prior),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_silent(waic(ao_bsir_with_prior)) 

# (4) method  = "HMC"   
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed))  
ao_hmc <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)
ao_hmc_two_chains <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 2L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)                       
expect_silent(ao_hmc)
expect_inherits(ao_hmc,"remstimate")  
expect_length(ao_hmc,2)
expect_identical(names(ao_hmc),c("sender_model","receiver_model"))
expect_length(ao_hmc$sender_model,10)
expect_length(ao_hmc$receiver_model,10)
expect_identical(names(ao_hmc$sender_model),c("draws","log_posterior","coefficients","post.mean","vcov","sd","loglik","df.null","df.model","df.residual"))
expect_identical(names(ao_hmc$receiver_model),c("draws","log_posterior","coefficients","post.mean","vcov","sd","loglik","df.null","df.model","df.residual"))
expect_length(attributes(ao_hmc),16)
expect_identical(names(attributes(ao_hmc)),c("names","class","formula","model","ordinal","method","approach","statistics","nsim","seed","nchains","burnin","thin","init","where_is_baseline","ncores"))
expect_identical(attr(ao_hmc,"approach"),"Bayesian")
expect_silent(print(ao_hmc))
expect_silent(summary(ao_hmc))
expect_silent(diagnostics(object = ao_hmc, reh = ao_reh, stats = ao_reh_stats))
ao_reh_diagnostics <- diagnostics(object = ao_hmc, reh = ao_reh, stats = ao_reh_stats)
expect_silent(plot(x = ao_hmc,reh = ao_reh, diagnostics = ao_reh_diagnostics))
ao_reh_diagnostics <- diagnostics(object = ao_hmc_two_chains, reh = ao_reh, stats = ao_reh_stats) # two chains
expect_silent(plot(x = ao_hmc_two_chains,reh = ao_reh, diagnostics = ao_reh_diagnostics)) # two chains
expect_error(aic(ao_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(aicc(ao_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_error(bic(ao_hmc),
"'approach' must be 'Frequentist'",
fixed = TRUE)
expect_silent(waic(ao_hmc)) 

# tests on "WAIC" for "HMC"                      
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 100L,
                        burnin = 5L,
                        seed = set_seed,
                        WAIC = TRUE,
                        nsimWAIC = 10))

# nsim = NULL
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = NULL,
                        burnin = 5L,
                        seed = set_seed))  

# L = NULL
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        L = NULL,
                        seed = set_seed))  

# epsilon = NULL
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        epsilon = NULL,
                        seed = set_seed))  

# ordinal likelihood (actor-oriented modeling)
ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor", ordinal = TRUE)
ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model, method="pt")
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "MLE"))


# Risk set "active"

# testing estimation methods with active riskset 
ao_reh <- remify::remify(edgelist = ao_data$edgelist, model = "actor", riskset = "active")

# calculating statistics
ao_reh_stats <- remstats::remstats(reh = ao_reh, sender_effects = rate_model, receiver_effects = choice_model, method="pt")
attr(ao_reh_stats,"formula") <- NULL

# (1) method = "MLE"
rate_model_two_effects <- ~ 1 + remstats::indegreeSender() + remstats::outdegreeSender()
ao_reh_stats_2 <- remstats::remstats(reh = ao_reh, sender_effects = rate_model_two_effects, receiver_effects = choice_model, method="pt")
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats_2,
                        ncores = 1L,
                        method = "MLE"))

ao_mle <- remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats_2,
                        ncores = 1L,
                        method = "MLE")
expect_silent(diagnostics(object = ao_mle, reh = ao_reh, stats = ao_reh_stats_2))

# (2) method = "GDADAMAX" 
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "GDADAMAX",
                        epochs = 10L))
# (3) method  = "BSIR" 

## (3.1) with prior = NULL
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = NULL,
                        seed = set_seed))     

## (3.2) with a specified prior 
priormvt <- mvnfast::dmvt
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "BSIR",
                        nsim = 10L,
                        prior = list(sender_model = priormvt, receiver_model = priormvt),
                        prior_args = list(sender_model =  list(mu = rep(0,dim(ao_reh_stats$sender_stats)[3]),
                                                            sigma = diag(dim(ao_reh_stats$sender_stats)[3]),
                                                            df = 1), 
                                        receiver_model = list(mu = rep(0,dim(ao_reh_stats$receiver_stats)[3]),
                                                            sigma = diag(dim(ao_reh_stats$receiver_stats)[3]),
                                                            df = 1)),
                        log = TRUE,
                        seed = set_seed
                        ))
# (4) method  = "HMC"    
expect_silent(remstimate::remstimate(reh = ao_reh,
                        stats = ao_reh_stats,
                        ncores = 1L,
                        method = "HMC",
                        nchains = 1L,
                        nsim = 10L,
                        burnin = 5L,
                        seed = set_seed)) 




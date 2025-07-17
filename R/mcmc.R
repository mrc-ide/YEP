# TODO: Replace the contents of this file with new functions for running fitting using monty (in new file)
# Global variable - additional parameter names
extra_param_names <- c("vaccine_efficacy","p_severe_inf","p_death_severe_inf","p_rep_severe",
                       "p_rep_death","m_FOI_Brazil")
#TODO - Add provision for parameters of temperature/precipitation functional forms?
#-------------------------------------------------------------------------------
#' @title MCMC
#'
#' @description Combined MCMC Multi-Region - series of MCMC iterations for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters.
#'
#' @param params_data #Data frame of parameter information containing names, initial values, maximum and minimum values,
#'   mean and standard deviation (for prior calculation) and flag indicating whether parameter estimated or fixed\cr
#'   Parameters to include: coefficients of environmental covariates to calculate FOI_spillover and R0, reported
#'   vaccination effectiveness, probability of severe case reporting, probability of fatal case reporting, Brazil
#'   FOI_spillover multiplier, FOI_spillover and R0 (latter two never estimated, included only for priors)\cr
#'   TBA - instructions
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives (TBA - instructions)
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths (TBA - instructions)
#' @param filename_prefix Prefix of output RDS file name, e.g. "Chain.Rds"
#' @param Niter Total number of iterations to run
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#' @param time_inc time increment in days (must be 1 or 5)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each iteration
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing time-varying environmental covariate data:\cr
#'   regions: Vector of region labels\cr
#'   env_vars: Vector of covariate names\cr
#'   values: Array of covariate values with dimensions (number of covariates, number of regions, number of time points).
#'   Number of time points must be correct for mode_time setting.\cr
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'   If mode_time = 0, no time variation (constant values)\cr
#'   If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'   If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'   If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'   If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'   If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' '
#' @export
#'
MCMC <- function(params_data = data.frame(name="FOI_var1",initial=1,max=Inf,min=-Inf,mean=0,sd=1,estimate=TRUE),
                 input_data = list(), obs_sero_data = NULL, obs_case_data = NULL, filename_prefix = "Chain", Niter = 1,
                 mode_start = 1, time_inc = 1.0, n_reps = 1, enviro_data_const = list(), enviro_data_var = list(),
                 deterministic = FALSE, mode_time = 1, mode_parallel = FALSE, cluster = NULL){


  #Section identical to equivalent section in mcmc_prelim_fit
  {
    # Checks
    assert_that(is.logical(deterministic))
    assert_that(mode_start %in% c(0, 1), msg = "mode_start must have value 0 or 1")
    if(is.null(obs_case_data)==FALSE){
      assert_that(all(obs_case_data$cases==round(obs_case_data$cases,0)),msg="Case data values must be integers")
      assert_that(all(obs_case_data$deaths==round(obs_case_data$deaths,0)),msg="Case data values must be integers")
    }
    regions = regions_breakdown(c(obs_sero_data$region,obs_case_data$region))
    assert_that(all(regions %in% enviro_data_const$region),
                msg = "Time-invariant environmental data must be available for all regions in observed data")
    if(is.null(enviro_data_var)==FALSE){
      assert_that(enviro_data_var_check(enviro_data_var))
      assert_that(all(regions %in% enviro_data_var$regions),
                  msg = "Time-variant environmental data must be available for all regions in observed data")
    }

    #Truncate input and environmental data to only include relevant regions
    input_data = input_data_truncate(input_data,regions)
    enviro_data_const = subset(enviro_data_const, enviro_data_const$region %in% regions)
    if(is.null(enviro_data_var)==FALSE){enviro_data_var = enviro_data_var_truncate(enviro_data_var,regions)}

    #Designate constant and variable covariates and check parameter data
    const_covars = colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
    var_covars = enviro_data_var$env_vars
    covar_names = c(const_covars,var_covars)
    checks <- mcmc_checks(params_data, covar_names, check_initial = TRUE)
    n_env_vars = length(covar_names)
    i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars]
    i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars]
    i_R0_const = i_FOI_const + n_env_vars
    i_R0_var = i_FOI_var + n_env_vars
    i_FOI_prior=c(1:nrow(params_data))[params_data$name=="FOI"]
    i_R0_prior=c(1:nrow(params_data))[params_data$name=="R0"]

    # Get names of extra estimated parameters (estimated parameters other than coefficients of environmental covariates)
    all_est_param_names = params_data$name[params_data$estimate==TRUE]
    n_params=length(all_est_param_names)
    add_est_param_names=all_est_param_names[c(((2*n_env_vars)+1):n_params)]

    # Cross-reference regions
    if(is.null(obs_sero_data)){xref_sero=NULL}else{xref_sero=template_region_xref(obs_sero_data,input_data$region_labels)}
    if(is.null(obs_case_data)){xref_case=NULL}else{xref_case=template_region_xref(obs_case_data,input_data$region_labels)}
  }

  #MCMC setup
  chain = chain_prop = posterior_current = posterior_prop = flag_accept = chain_cov_all = NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  log_params = log(params_data$initial[params_data$estimate==TRUE])
  names(log_params) = all_est_param_names
  chain_cov = 1
  adapt = 0
  posterior_value_current = -Inf

  #time0=Sys.time()
  #Iterative estimation
  for (iter in 1:Niter){

    #Propose new parameter values
    log_params_prop = param_prop_setup(log_params, chain_cov, adapt)

    #Calculate likelihood using single_posterior_calc function
    posterior_value_prop = single_posterior_calc(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                                 params_data = params_data, mode_start = mode_start, time_inc = time_inc,
                                                 n_reps = n_reps, enviro_data_const = enviro_data_const,
                                                 enviro_data_var = enviro_data_var, add_est_param_names = add_est_param_names,
                                                 deterministic = deterministic, mode_time = mode_time,
                                                 mode_parallel = mode_parallel, cluster = cluster,
                                                 i_FOI_const = i_FOI_const, i_FOI_var = i_FOI_var,
                                                 i_R0_const = i_R0_const, i_R0_var = i_R0_var,
                                                 i_FOI_prior = i_FOI_prior, i_R0_prior = i_R0_prior,
                                                 n_env_vars = n_env_vars, xref_sero = xref_sero, xref_case = xref_case)
    gc() #Clear garbage to prevent memory creep

    if(is.finite(posterior_value_prop) == FALSE) {
      p_accept = -Inf
    } else {
      p_accept = posterior_value_prop - posterior_value_current
      if(is.na(p_accept) ){ p_accept = -Inf}
    }

    ## accept/reject iteration:
    tmp = runif(1)
    if(tmp<min(exp(p_accept), 1)) {
      log_params = log_params_prop
      posterior_value_current = posterior_value_prop
      accept = 1
    } else {accept = 0}

    #save current iteration
    chain = rbind(chain, log_params)
    chain_prop = rbind(chain_prop, log_params_prop)
    posterior_current = rbind(posterior_current, posterior_value_current)
    posterior_prop = rbind(posterior_prop, posterior_value_prop)
    flag_accept = rbind(flag_accept, accept)
    chain_cov_all = rbind(chain_cov_all, max(chain_cov))

    #Set output headings
    if(iter == 1){
      colnames(chain) = colnames(chain_prop) = all_est_param_names
      for(i in 1:n_params){colnames(chain_prop)[i] = paste("Test_", colnames(chain_prop)[i], sep = "")}
      colnames(posterior_current) = "posterior_current"
      colnames(posterior_prop) = "posterior_prop"
      colnames(flag_accept) = "flag_accept"
      colnames(chain_cov_all) = "chain_cov_all"
    }

    #Output chain to file every 10 iterations
    if (iter %% 10 == 0){
      # if (iter %% 10000 == 10){fileIndex = (iter-10)/10000}
      # if(fileIndex >=  10){fn = paste0(filename_prefix,fileIndex,".csv")}else{fn = paste0(filename_prefix,"0",fileIndex,".csv")}
      fn=paste0(filename_prefix,".Rds")
      if(file.exists(fn) == FALSE){file.create(fn)}
      # lines = min(((fileIndex*10000) + 1), iter):iter

      try(saveRDS(data.frame(cbind(posterior_current, posterior_prop, exp(chain), flag_accept,
                                   exp(chain_prop), chain_cov_all)[c(1:iter),]),file = fn))
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov = cov(chain[max(nrow(chain)-10000, 1):nrow(chain), ])
    } else {
      adapt = 0
      chain_cov = 1
    }
    # time1=Sys.time()
    # cat("\n Iteration ",iter,"\tTime elapsed:",time1-time0,"\tRate: ",iter/(as.numeric(time1)-as.numeric(time0)))
  }

  return(data.frame(cbind(posterior_current, exp(chain))))
  #return(NULL)
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit
#'
#' @description Test multiple sets of parameters randomly drawn from range between maximum and minimum
#' values in order to find approximate values giving maximum posterior likelihood
#'
#' @details This function is used to estimate the model parameter values giving maximum posterior likelihood; it is
#' primarily intended to be used to generate initial parameter values for Markov Chain Monte Carlo fitting (using
#' the mcmc() function).
#'
#' @param n_iterations = Number of times to run and adjust maximum/minimum
#' @param n_param_sets = Number of parameter sets to run in each iteration
#' @param n_bounds = Number of parameter sets (with highest likelihood values) to take at each iteration to create new
#' maximum/minimum values
#' @param params_data TBA
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#' @param time_inc time increment in days (must be 1 or 5)
#' @param n_reps Number of repetitions
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing values of time-varying environmental covariates (TBA)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_time TBA
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' @param plot_graphs TRUE/FALSE - plot graphs of evolving parameter space
#' '
#' @export
#'
mcmc_prelim_fit <- function(n_iterations = 1, n_param_sets = 1, n_bounds = 1, params_data = list(),
                            input_data = list(), obs_sero_data = list(), obs_case_data = list(),
                            mode_start = 1, time_inc = 1.0, n_reps = 1, enviro_data_const = list(), enviro_data_var = list(),
                            deterministic = TRUE, mode_time = 0, mode_parallel = FALSE, cluster = NULL, plot_graphs = FALSE){

  #Section identical to equivalent section in MCMC
  {
    # Checks
    assert_that(is.logical(deterministic))
    assert_that(mode_start %in% c(0, 1), msg = "mode_start must have value 0 or 1")
    if(is.null(obs_case_data)==FALSE){
      assert_that(all(obs_case_data$cases==round(obs_case_data$cases,0)),msg="Case data values must be integers")
      assert_that(all(obs_case_data$deaths==round(obs_case_data$deaths,0)),msg="Case data values must be integers")
    }
    regions = regions_breakdown(c(obs_sero_data$region,obs_case_data$region))
    assert_that(all(regions %in% enviro_data_const$region),
                msg = "Time-invariant environmental data must be available for all regions in observed data")
    if(is.null(enviro_data_var)==FALSE){
      assert_that(enviro_data_var_check(enviro_data_var))
      assert_that(all(regions %in% enviro_data_var$regions),
                  msg = "Time-variant environmental data must be available for all regions in observed data")
    }

    #Truncate input and environmental data to only include relevant regions
    input_data = input_data_truncate(input_data,regions)
    enviro_data_const = subset(enviro_data_const, enviro_data_const$region %in% regions)
    if(is.null(enviro_data_var)==FALSE){enviro_data_var = enviro_data_var_truncate(enviro_data_var,regions)}

    #Designate constant and variable covariates and check parameter data
    const_covars = colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
    var_covars = enviro_data_var$env_vars
    covar_names = c(const_covars,var_covars)
    checks <- mcmc_checks(params_data, covar_names, check_initial = FALSE)
    n_env_vars = length(covar_names)
    i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars]
    i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars]
    i_R0_const = i_FOI_const + n_env_vars
    i_R0_var = i_FOI_var + n_env_vars
    i_FOI_prior=c(1:nrow(params_data))[params_data$name=="FOI"]
    i_R0_prior=c(1:nrow(params_data))[params_data$name=="R0"]

    # Get names of extra estimated parameters (estimated parameters other than coefficients of environmental covariates)
    all_est_param_names = params_data$name[params_data$estimate==TRUE]
    n_params=length(all_est_param_names)
    add_est_param_names=all_est_param_names[c(((2*n_env_vars)+1):n_params)]

    # Cross-reference regions
    if(is.null(obs_sero_data)){xref_sero=NULL}else{xref_sero=template_region_xref(obs_sero_data,input_data$region_labels)}
    if(is.null(obs_case_data)){xref_case=NULL}else{xref_case=template_region_xref(obs_case_data,input_data$region_labels)}
  }

  best_fit_results = list()
  log_params_min=log(params_data$min[params_data$estimate==TRUE])
  log_params_max=log(params_data$max[params_data$estimate==TRUE])
  assert_that(all(is.finite(c(log_params_min,log_params_max))),msg="Limits must be finite")

  if(plot_graphs){
    xlabels = all_est_param_names
    for(i in 1:n_params){xlabels[i] = substr(xlabels[i], 1, 15)}
    ylabels = 10^c(-8, -6, -4, -3, -2, -1, 0, 1)
    par(mar = c(6, 2, 1, 1))
    ylim = c(min(log_params_min), max(log_params_max))
  }

  for(iteration in 1:n_iterations){
    cat("\nIteration: ", iteration, "\n", sep = "")
    all_param_sets <- lhs(n = n_param_sets, rect = cbind(log_params_min, log_params_max))
    results = data.frame()

    for(set in 1:n_param_sets){
      cat("\n\tSet: ", set, sep = "")
      log_params_prop = all_param_sets[set, ]

      cat("\n\tParams: ", signif(log_params_prop, 3))

      names(log_params_prop) = all_est_param_names
      posterior_value = single_posterior_calc(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                              params_data = params_data, mode_start = mode_start, time_inc = time_inc,
                                              n_reps = n_reps, enviro_data_const = enviro_data_const,
                                              enviro_data_var = enviro_data_var, add_est_param_names = add_est_param_names,
                                              deterministic = deterministic, mode_time = mode_time,
                                              mode_parallel = mode_parallel, cluster = cluster,
                                              i_FOI_const = i_FOI_const, i_FOI_var = i_FOI_var,
                                              i_R0_const = i_R0_const, i_R0_var = i_R0_var,
                                              i_FOI_prior = i_FOI_prior, i_R0_prior = i_R0_prior,
                                              n_env_vars = n_env_vars, xref_sero = xref_sero, xref_case = xref_case)
      gc() #Clear garbage to prevent memory creep
      results <- rbind(results, c(set, exp(log_params_prop), posterior_value))
      if(set == 1){colnames(results) = c("set", all_est_param_names, "posterior")}

      cat("\n\tPosterior likelihood = ", posterior_value, sep = "")

    }
    results <- results[order(results$posterior, decreasing = TRUE), ]
    best_fit_results[[iteration]] = results

    log_params_min_new = log_params_max_new = rep(0, n_params)
    for(i in 1:n_params){
      log_params_min_new[i] = min(log(results[c(1:n_bounds), i + 1]))
      log_params_max_new[i] = max(log(results[c(1:n_bounds), i + 1]))
    }
    names(log_params_min_new) = names(log_params_max_new) = all_est_param_names

    if(plot_graphs){
      matplot(x = c(1:n_params), y = log(t(results[c(1:n_bounds), c(1:n_params) + 1])), type = "p", pch = 16, col = 1,
              xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = ylim)
      axis(side = 1, at = c(1:n_params), labels = xlabels, las = 2, cex.axis = 0.7)
      axis(side = 2, at = log(ylabels), labels = ylabels)
      matplot(x = c(1:n_params), y = log_params_min, type = "l", col = 1, lty = 2, add = TRUE)
      matplot(x = c(1:n_params), y = log_params_max, type = "l", col = 1, lty = 2, add = TRUE)
      matplot(x = c(1:n_params), y = log_params_min_new, type = "l", col = 2, add = TRUE)
      matplot(x = c(1:n_params), y = log_params_max_new, type = "l", col = 2, add = TRUE)
    }

    log_params_min = log_params_min_new
    log_params_max = log_params_max_new
  }

  return(best_fit_results)
}
#-------------------------------------------------------------------------------
#' @title single_posterior_calc
#'
#' @description Function which calculates and outputs posterior likelihood of observing simulated data
#'
#' @details This function calculates the posterior likelihood of observing a set of observations (across multiple
#' regions and data types) for a given proposed parameter set. [TBA]
#'
#' @param log_params_prop Proposed values of parameters to be estimated (natural logarithm of actual parameters)
#' @param input_data List of population and vaccination data for multiple regions (created using data input
#'   creation code and usually loaded from RDS file), with cross-reference tables added in MCMC
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param ... = Constant parameters/flags/etc. loaded to or determined by mcmc() and mcmc_prelim_fit, including params_data,
#'   mode_start, time_inc, n_reps, enviro_data, add_est_param_names,
#'   deterministic, mode_time, mode_parallel, cluster, i_FOI_const, i_FOI_var, i_R0_const, i_R0_var, i_FOI_prior, i_R0_prior
#'
#' @export
#'
single_posterior_calc <- function(log_params_prop = c(),input_data = list(),obs_sero_data = NULL,obs_case_data = NULL,...){

  consts = list(...)

  #Get additional values, calculate associated prior values
  prior_like = 0
  vaccine_efficacy = p_severe_inf = p_death_severe_inf = p_rep_severe = p_rep_death = m_FOI_Brazil = 1.0
  #TODO - Add provision for parameters of temperature/precipitation functional forms?
  for(var_name in extra_param_names){
    i=match(var_name,consts$params_data$name)
    if(var_name %in% consts$add_est_param_names){
      value = exp(as.numeric(log_params_prop[match(var_name, names(log_params_prop))]))
      assign(var_name, value)
      prior_like = prior_like + log(dtrunc(value, "norm", a = consts$params_data$min[i], b = consts$params_data$max[i],
                                           mean = consts$params_data$mean[i], sd = consts$params_data$sd[i]))
    } else {
      assign(var_name, consts$params_data$initial[i])
    }
  }

  #If prior is finite so far, get normal-distribution prior values for environmental coefficients
  if(is.finite(prior_like)){
    for(i in 1:(2*consts$n_env_vars)){
      prior_like = prior_like + log(dtrunc(log_params_prop[i], "norm",
                                           a = log(consts$params_data$min[i]), b = log(consts$params_data$max[i]),
                                           mean = log(consts$params_data$mean[i]), sd = log(consts$params_data$sd[i])))
    }
  }

  #If prior is finite so far, get FOI and R0 values and calculate any associated prior
  if(is.finite(prior_like)){

    #TODO - If parameters of temperature/precipitation functional forms included, add functional form conversion of temp/precip enviro data?
    #consts$enviro_data_var

    FOI_values = epi_param_calc(exp(log_params_prop[consts$i_FOI_const]), exp(log_params_prop[consts$i_FOI_var]),
                                consts$enviro_data_const, consts$enviro_data_var)
    R0_values = epi_param_calc(exp(log_params_prop[consts$i_R0_const]), exp(log_params_prop[consts$i_R0_var]),
                               consts$enviro_data_const, consts$enviro_data_var)

    for(n_region in 1:length(input_data$region_labels)){ #Apply Brazil FOI multiplier to relevant regions
      if(substr(input_data$region_labels[n_region],1,3) == "BRA"){FOI_values[n_region] = FOI_values[n_region]*m_FOI_Brazil}
    }
    prior_like = prior_like  +
      sum(log(dtrunc(rowMeans(FOI_values), "norm", #TODO - change to use log(FOI)?
                     a = consts$params_data$min[consts$i_FOI_prior], b = consts$params_data$max[consts$i_FOI_prior],
                     mean = consts$params_data$mean[consts$i_FOI_prior], sd = consts$params_data$sd[consts$i_FOI_prior]))) +
      sum(log(dtrunc(rowMeans(R0_values), "norm",
                     a = consts$params_data$min[consts$i_R0_prior], b = consts$params_data$max[consts$i_R0_prior],
                     mean = consts$params_data$mean[consts$i_R0_prior], sd = consts$params_data$sd[consts$i_R0_prior])))
  }

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_like)) {

    #Generate modelled data over all regions
    dataset <- Generate_Dataset(FOI_values, R0_values, input_data, obs_sero_data, obs_case_data, vaccine_efficacy,
                                consts$time_inc,consts$mode_start, start_SEIRV = NULL, consts$mode_time, consts$n_reps,
                                consts$deterministic, p_severe_inf, p_death_severe_inf, p_rep_severe, p_rep_death,
                                consts$mode_parallel, consts$cluster, output_frame = FALSE, seed = NULL, #TBC?
                                xref_sero = consts$xref_sero, xref_case = consts$xref_case)

    #Likelihood of observing serological data
    if(is.null(obs_sero_data) == FALSE){
      sero_like_values = sero_data_compare(dataset$model_sero_values, obs_sero_data)
    } else {sero_like_values = 0}

    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data) == FALSE){
      cases_like_values = case_data_compare(dataset$model_case_values, obs_case_data$cases)
      if(is.null(obs_case_data$deaths) == FALSE){
        deaths_like_values = case_data_compare(dataset$model_death_values, obs_case_data$deaths)
      } else {deaths_like_values = 0}
    } else {cases_like_values = deaths_like_values = 0}

    posterior = prior_like + sum(sero_like_values, na.rm = TRUE) + sum(cases_like_values, na.rm = TRUE) +
      sum(deaths_like_values, na.rm = TRUE)

  } else {posterior = -Inf}

  return(posterior)
}
#-------------------------------------------------------------------------------
#' @title mcmc_checks
#'
#' @description Perform checks on parameter data
#'
#' @details This function, which is called by MCMC() and mcmc_prelim_fit(), performs a number of checks on parameter data
#' to ensure proper functionality. [TBA]
#'
#' @param params_data TBA
#' @param covar_names TBA
#' @param check_initial TBA
#'
#' @export
#'
mcmc_checks <- function(params_data = list(), covar_names = c(), check_initial = FALSE){

  #Check columns
  assert_that(all(colnames(params_data)==c("name","initial","max","min","mean","sd","estimate")))

  # Check number of parameters
  n_env_vars = length(covar_names)
  n_rows = nrow(params_data)
  assert_that(n_rows==(2*n_env_vars) + length(extra_param_names) + 2, msg = "Check number of rows of parameter data")
  range1=c(1:(n_rows-2))
  range2=c(((2*n_env_vars)+1):(n_rows-2))

  # Check parameter names
  assert_that(all(params_data$name[c(1:n_env_vars)] == paste0("FOI_",covar_names)),
              msg = "Check covariate FOI coefficient parameter names")
  assert_that(all(params_data$name[c(1:n_env_vars) + n_env_vars] == paste0("R0_", covar_names)),
              msg = "Check covariate R0 coefficient parameter names")
  assert_that(all(params_data$name[c(((2*n_env_vars)+1):n_rows)]==c(extra_param_names,"FOI","R0")),
              msg = "Check additional parameter names")

  # Check estimation flags
  assert_that(all(params_data$estimate[c(1:(2*n_env_vars))]), msg = "All covariate FOI/R0 coefficients must be estimated")
  assert_that(all(params_data$estimate[c((n_rows-1):n_rows)]==FALSE), msg = "FOI/R0 values never estimated, priors only")

  # Check parameter max-min ranges and mean values
  assert_that(all(params_data$max>=params_data$min), msg = "Check parameter max-min ranges")
  assert_that(all(params_data$mean>=params_data$min), msg = "Check parameter mean and minimum values")
  assert_that(all(params_data$mean<=params_data$max), msg = "Check parameter mean and maximum values")
  assert_that(all(params_data$min[range2]>=0), msg = "Check additional parameter minimum values")
  assert_that(all(params_data$max[range2]<=1), msg = "Check additional parameter maximum values")

  # If relevant, check initial values within max-min ranges
  if(check_initial){
    assert_that(all(params_data$initial[range1]>=params_data$min[range1]), msg = "Check parameter initial and minimum values")
    assert_that(all(params_data$initial[range1]<=params_data$max[range1]), msg = "Check parameter initial and maximum values")
  }

  # Check parameter standard deviations
  assert_that(all(params_data$sd[params_data$estimate==TRUE]>0), msg = "Check parameter standard deviations")

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title mcmc_params_data_create
#'
#' @description Create data frame of parameter data for input into MCMC() and mcmc_prelim_fit()
#'
#' @details TBA
#'
#' @param covar_names TBA
#' @param add_est_param_names TBA
#'
#' @export
#'
mcmc_params_data_create <- function(covar_names = c("Var1", "Var2"),
                                    add_est_param_names = "vaccine_efficacy"){
  assert_that(all(c("vaccine_efficacy","p_severe_inf","p_death_severe_inf","p_rep_severe",
                    "p_rep_death") %in% extra_param_names))
  assert_that(all(add_est_param_names %in% extra_param_names))

  n_env_vars=length(covar_names)
  n_extra=length(extra_param_names)
  nrows=(2*n_env_vars)+n_extra+2

  #Create data frame
  params_data=data.frame(name=c(paste0("FOI_",covar_names),paste0("R0_",covar_names),extra_param_names,
                                "FOI","R0"),
                         initial=c(rep(1e-5,n_env_vars),rep(1e-3,n_env_vars),rep(1,n_extra),NA,NA),
                         max=c(rep(1e-4,n_env_vars),rep(10,n_env_vars),rep(1,n_extra),Inf,Inf),
                         min=c(rep(1e-10,n_env_vars),rep(1e-4,n_env_vars),rep(1e-3,n_extra),0,0),
                         mean=c(rep(1e-10,n_env_vars),rep(1,n_env_vars),rep(1,n_extra),0,4.0),
                         sd=rep(1000,nrows),
                         estimate=c(rep(TRUE,2*n_env_vars),rep(FALSE,n_extra+2)))

  #Set additional parameters to be estimated
  params_data$estimate[params_data$name %in% add_est_param_names] = TRUE

  #Set values for vaccine_efficacy
  # i_ve=c(1:nrows)[params_data$name=="vaccine_efficacy"]
  # params_data$initial

  #Set values for p_severe_inf
  i_psi=c(1:nrows)[params_data$name=="p_severe_inf"]
  params_data$initial[i_psi]=params_data$max[i_psi]=params_data$min[i_psi]=params_data$mean[i_psi]=0.12

  #Set values for p_death_severe_inf
  i_pdsi=c(1:nrows)[params_data$name=="p_death_severe_inf"]
  params_data$initial[i_pdsi]=params_data$max[i_pdsi]=params_data$min[i_pdsi]=params_data$mean[i_pdsi]=0.39

  return(params_data)
}
#-------------------------------------------------------------------------------
#' @title param_prop_setup
#'
#' @description Set up proposed new log parameter values for next iteration in chain
#'
#' @details Takes in current values of parameter set used for Markov Chain Monte Carlo fitting and proposes new values
#' from multivariate normal distribution where the existing values form the mean and the standard deviation is
#' based on the chain covariance or (if the flag "adapt" is set to 1) a flat value based on the number of parameters.
#'
#' @param log_params Previous log parameter values used as input
#' @param chain_cov Covariance calculated from previous iterations in chain
#' @param adapt 0/1 flag indicating which type of covariance to use for proposition value (TBA)
#' '
#' @export
#'
param_prop_setup <- function(log_params = c(), chain_cov = 1, adapt = 0){

  n_params = length(log_params)
  if (adapt == 1) {
    sigma = (5.6644*chain_cov)/n_params #'optimal' scaling of chain covariance (2.38 ^ 2)
    log_params_prop_a = rmvnorm(n = 1, mean = log_params, sigma = sigma)
  } else {
    sigma = (1.0e-4*diag(n_params))/n_params #this is an inital proposal covariance, see [Mckinley et al 2014] ((1e-2) ^ 2)
    log_params_prop_a = rmvnorm(n = 1, mean = log_params, sigma = sigma)
  }

  return(log_params_prop_a[1, ])
}

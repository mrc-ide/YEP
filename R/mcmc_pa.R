# Version of mcmc.R with functions updated to streamline parameter implementation
#-------------------------------------------------------------------------------
#' @title MCMC3
#'
#' @description Combined MCMC Multi-Region - series of MCMC iterations for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters.
#'
#' @param params_data #Data frame of parameter information containing names, initial values and optionally:\cr
#'   Maximum and minimum values (for non-zero prior)\cr
#'   Mean and standard deviation (for normal prior)
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param filename_prefix Prefix of names for output files; function outputs a CSV file every 10,000 iterations with a
#'   name in the format: "(filename_prefix)XX.csv", e.g. Chain00.csv
#' @param Niter Total number of iterations to run
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param prior_settings List containing settings for priors: must contain text named "type":
#'   If type = "zero", prior probability is always zero \cr
#'   If type = "flat", prior probability is zero if log parameter values in designated maximum/minimum range
#'   given in params_data\cr
#'   If type = "norm", prior probability is given by truncated normal distribution using limits, mean and standard
#'   deviation values in params_data and additional values: \cr
#'   + FOI_mean + FOI_sd (mean + standard deviation of computed FOI, single values)  \cr
#'   + R0_mean + R0_sd (mean + standard deviation of computed R0, single values) \cr
#' @param time_inc time increment in days (must be 1 or 5)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each iteration
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing time-varying environmental covariate data:\cr
#'   regions: Vector of region labels\cr
#'   env_vars: Vector of covariate names\cr
#'   values: Array of covariate values with dimensions (number of covariates, number of regions, number of time points).
#'   Number of time points must be correct for mode_time setting.\cr
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set \cr
#'   vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'   p_rep_severe Probability of observation of severe infection \cr
#'   p_rep_death Probability of observation of death \cr
#'   m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
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
MCMC3 <- function(params_data = c(), input_data = list(), obs_sero_data = NULL, obs_case_data = NULL,
                 filename_prefix = "Chain", Niter = 1, mode_start = 1, prior_settings = list(type = "zero"),
                 time_inc = 1.0, n_reps = 1, enviro_data_const = list(), enviro_data_var = list(),
                 p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                 add_values = list(vaccine_efficacy = 1.0, p_rep_severe = 1.0, p_rep_death = 1.0, m_FOI_Brazil = 1.0),
                 deterministic = FALSE, mode_time = 1, mode_parallel = FALSE, cluster = NULL){

  assert_that(is.logical(deterministic))
  assert_that(mode_start %in% c(0, 1, 3), msg = "mode_start must have value 0, 1 or 3 (NB 3 should be changed to 1)")
  if(is.null(obs_case_data) == FALSE){assert_that(all(c("p_rep_severe","p_rep_death") %in% names(add_values)),
                                                  msg = "Reporting probabilities required for case data")}

  #Process input data to check that all regions with sero and/or case data supplied are present, remove
  #regions without any supplied data.
  regions = regions_breakdown(c(obs_sero_data$region,obs_case_data$region))
  input_data = input_data_truncate(input_data,regions)
  n_regions = length(input_data$region_labels)

  # Take subset of environmental data and check that environmental data available for all regions
  assert_that(all(regions %in% enviro_data_const$region),
              msg = "Time-invariant environmental data must be available for all regions in observed data")
  enviro_data_const = subset(enviro_data_const, enviro_data_const$region %in% regions)
  if(is.null(enviro_data_var)==FALSE){
    assert_that(enviro_data_var_check(enviro_data_var))
    assert_that(all(regions %in% enviro_data_var$regions),
                msg = "Time-variant environmental data must be available for all regions in observed data")
    enviro_data_var = enviro_data_var_truncate(enviro_data_var,regions)
  }

  #Get names of additional parameters to be estimated
  extra_estimated_params = c()
  for(var_name in names(add_values)){
    if(is.na(add_values[[var_name]]) == TRUE){extra_estimated_params = append(extra_estimated_params, var_name)}
  }

  #Run checks on inputs
  checks <- mcmc_checks3(params_data, n_regions, prior_settings, enviro_data_const, enviro_data_var, add_values,
                        extra_estimated_params)

  #Designate constant and variable covariates
  const_covars = colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
  var_covars = enviro_data_var$env_vars
  covar_names = c(const_covars,var_covars)
  n_env_vars = length(covar_names)
  i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars]
  i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars]
  i_R0_const = i_FOI_const + n_env_vars
  i_R0_var = i_FOI_var + n_env_vars

  #MCMC setup
  chain = chain_prop = posterior_current = posterior_prop = flag_accept = chain_cov_all = NULL
  n_params = nrow(params_data)
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  log_params = params_data$initial
  chain_cov = 1
  adapt = 0
  posterior_value_current = -Inf

  #Iterative estimation
  for (iter in 1:Niter){

    #Propose new parameter values
    log_params_prop = param_prop_setup(log_params, chain_cov, adapt)

    #Calculate likelihood using single_posterior_calc3 function
    posterior_value_prop = single_posterior_calc3(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                                 params_data = params_data, mode_start = mode_start,prior_settings = prior_settings,
                                                 time_inc = time_inc, n_reps = n_reps, enviro_data_const = enviro_data_const,
                                                 enviro_data_var = enviro_data_var, p_severe_inf = p_severe_inf,
                                                 p_death_severe_inf = p_death_severe_inf, add_values = add_values,
                                                 extra_estimated_params = extra_estimated_params,
                                                 deterministic = deterministic, mode_time = mode_time,
                                                 mode_parallel = mode_parallel, cluster = cluster,
                                                 i_FOI_const = i_FOI_const, i_FOI_var = i_FOI_var,
                                                 i_R0_const = i_R0_const, i_R0_var = i_R0_var, n_env_vars = n_env_vars)
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
      colnames(chain) = colnames(chain_prop) = params_data$name
      for(i in 1:n_params){colnames(chain_prop)[i] = paste("Test_", colnames(chain_prop)[i], sep = "")}
      colnames(posterior_current) = "posterior_current"
      colnames(posterior_prop) = "posterior_prop"
      colnames(flag_accept) = "flag_accept"
      colnames(chain_cov_all) = "chain_cov_all"
    }

    #Output chain to file every 10 iterations; start new file every 10,000 iterations
    if (iter %% 10 == 0){
      if (iter %% 10000 == 10){fileIndex = (iter-10)/10000}
      if(fileIndex >=  10){fn = paste0(filename_prefix,fileIndex,".csv")}else{fn = paste0(filename_prefix,"0",fileIndex,".csv")}
      if(file.exists(fn) == FALSE){file.create(fn)}
      lines = min(((fileIndex*10000) + 1), iter):iter

      data_out <- cbind(posterior_current, posterior_prop, exp(chain), flag_accept, exp(chain_prop), chain_cov_all)[lines,]
      if(fileAccess(fn, 2) == 0){write.csv(data_out, fn, row.names = FALSE)}
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov = cov(chain[max(nrow(chain)-10000, 1):nrow(chain), ])
    } else {
      adapt = 0
      chain_cov = 1
    }
  }

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title single_posterior_calc3
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
#' @param ... = Constant parameters/flags/etc. loaded to or determined by mcmc() and mcmc_prelim_fit, including mode_start,
#' prior_settings, time_inc, n_reps, enviro_data, p_severe_inf, p_death_severe_inf, add_values, extra_estimated_params,
#' deterministic, mode_time, mode_parallel, cluster, i_FOI_const, i_FOI_var, i_R0_const, i_R0_var
#'
#' @export
#'
single_posterior_calc3 <- function(log_params_prop = c(),input_data = list(),obs_sero_data = NULL,obs_case_data = NULL,...){

  consts = list(...)

  #Check values for flat prior
  prior_like = 0
  if(consts$prior_settings$type == "flat"){
    if(any(log_params_prop<consts$params_data$min) ||
       any(log_params_prop>consts$params_data$max)){
      prior_like = -Inf}
  }

  #Get additional values, calculate associated normal-distribution prior values if relevant
  if(is.finite(prior_like)){
    vaccine_efficacy = p_rep_severe = p_rep_death = m_FOI_Brazil = 1.0
    for(var_name in names(consts$add_values)){
      if(var_name %in% consts$extra_estimated_params){
        i = match(var_name, names(log_params_prop))
        value = exp(as.numeric(log_params_prop[i]))
        assign(var_name, value)
        if(consts$prior_settings$type == "norm"){
          prior_like = prior_like + log(dtrunc(value, "norm", a = consts$params_data$min[i],
                                               b = consts$params_data$max[i],
                                               mean = consts$params_data$mean[i],
                                               sd = consts$params_data$sd[i]))
        }
      } else {assign(var_name, consts$add_values[[var_name]])}
    }
  }

  #If prior is finite so far, get normal-distribution prior values for environmental coefficients if relevant
  if(is.finite(prior_like) && consts$prior_settings$type == "norm"){
    for(i in 1:(2*consts$n_env_vars)){
      prior_like = prior_like + log(dtrunc(log_params_prop[i], "norm", a = log(consts$params_data$min[i]),
                                           b = log(consts$params_data$max[i]),
                                           mean = log(consts$params_data$mean[i]),
                                           sd = log(consts$params_data$sd[i])))
    }
  }

  #If prior is finite so far, get FOI and R0 values and calculate any associated prior
  if(is.finite(prior_like)){
    regions = input_data$region_labels
    n_regions = length(regions)

    FOI_values = epi_param_calc(coeffs_const = exp(log_params_prop[consts$i_FOI_const]),
                                coeffs_var = exp(log_params_prop[consts$i_FOI_var]),
                                enviro_data_const = consts$enviro_data_const,enviro_data_var = consts$enviro_data_var)
    R0_values = epi_param_calc(coeffs_const = exp(log_params_prop[consts$i_R0_const]),
                               coeffs_var = exp(log_params_prop[consts$i_R0_var]),
                               enviro_data_const = consts$enviro_data_const,enviro_data_var = consts$enviro_data_var)

    for(n_region in 1:n_regions){
      if(substr(regions[n_region],1,3) == "BRA"){FOI_values[n_region] = FOI_values[n_region]*m_FOI_Brazil}
    }
    if(consts$prior_settings$type == "norm"){
      prior_like = prior_like  +
        sum(log(dtrunc(rowMeans(R0_values),
                       "norm", a = 0, b = Inf, mean = consts$prior_settings$R0_mean, sd = consts$prior_settings$R0_sd))) +
        sum(log(dtrunc(rowMeans(FOI_values),
                       "norm", a = 0, b = 1, mean = consts$prior_settings$FOI_mean, sd = consts$prior_settings$FOI_sd)))
    }
  }

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_like)) {

    #Generate modelled data over all regions
    dataset <- Generate_Dataset(FOI_values, R0_values, input_data, obs_sero_data, obs_case_data, vaccine_efficacy,
                                consts$time_inc,consts$mode_start, start_SEIRV = NULL, consts$mode_time, consts$n_reps,
                                consts$deterministic, consts$p_severe_inf, consts$p_death_severe_inf, p_rep_severe, p_rep_death,
                                consts$mode_parallel, consts$cluster,output_frame = FALSE)

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
#' @title mcmc_checks3
#'
#' @description Perform checks on MCMC inputs
#'
#' @details This function, which is called by MCMC(), performs a number of checks on data to be used in fitting to
#' ensure proper functionality. It verifies that the number of parameters being estimated is consistent with other
#' settings and that certain values are not outwith sensible boundaries (e.g. probabilities must be between 0 and 1).
#'
#' @param params_data TBA
#' @param n_regions Number of regions
#' @param prior_settings List containing settings for priors; see documentation for MCMC() function for more details)
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing values of time-varying environmental covariates (TBA)
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the parameter set to be estimated \cr
#'  vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'  p_rep_severe Probability of observation of severe infection \cr
#'  p_rep_death Probability of observation of death \cr
#'  m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered) \cr
#' @param extra_estimated_params Vector of names of parameters to be estimated in addition to those governing FOI and R0;
#' see add_values
#'
#' @export
#'
mcmc_checks3 <- function(params_data = list(), n_regions = 1, prior_settings = list(type = "zero"),
                        enviro_data_const = list(), enviro_data_var = list(),
                        add_values = list(vaccine_efficacy = 1.0,p_rep_severe = 1.0,p_rep_death = 1.0,m_FOI_Brazil = 1.0),
                        extra_estimated_params = c()){

  n_params = nrow(params_data)
  assert_that(prior_settings$type %in% c("zero", "flat", "norm"),
              msg = "Prior settings type must be 'zero', 'flat' or 'norm'")
  if(prior_settings$type == "norm"){
    assert_that(is.numeric(prior_settings$R0_mean), msg = "Check prior_settings$R0_mean")
    assert_that(is.numeric(prior_settings$R0_sd), msg = "Check prior_settings$R0_sd")
    assert_that(is.numeric(prior_settings$FOI_mean), msg = "Check prior_settings$FOI_mean")
    assert_that(is.numeric(prior_settings$FOI_sd), msg = "Check prior_settings$FOI_sd")
  }

  # Check additional values
  add_value_names = names(add_values)
  assert_that("vaccine_efficacy" %in% add_value_names,
              msg = "Reported vaccination effectiveness vaccine_efficacy must be included in add_values")
  for(var_name in add_value_names){
    if(var_name %in% extra_estimated_params){
      assert_that(is.na(add_values[[var_name]]),msg = "Additional estimated parameters must be set to NA in add_values")
    } else {assert_that(add_values[[var_name]] <=  1.0 && add_values[[var_name]] >=  0.0,
                        msg = "Fixed additional parameters must be in range 0-1")}
  }

  # Get names of environmental covariates
  assert_that(is.null(enviro_data_const) == FALSE, msg = "Constant environmental data required")
  #assert_that(is.null(enviro_data_var) == FALSE, msg = "Variable environmental data required")
  covar_names = c(names(enviro_data_const[c(2:ncol(enviro_data_const))]),enviro_data_var$env_vars) #TBC
  n_env_vars = length(covar_names)

  # Check that total number of parameters is correct; check parameters named in correct order (TBA)
  if(is.null(extra_estimated_params)){n_extra_params = 0}else{n_extra_params = length(extra_estimated_params)}
  assert_that(n_params == (2*n_env_vars) + n_extra_params,
              msg = "Length of initial parameter vector must equal twice number of environmental covariates +
              number of additional estimated parameters")
  for(i in 1:n_env_vars){
    assert_that(params_data$name[i] == paste0("FOI_",covar_names[i]),msg = "Initial parameter vector must start with FOI coefficients")
    assert_that(params_data$name[i + n_env_vars] == paste0("R0_", covar_names[i]),
                msg = "R0 coefficients must follow FOI coefficients in initial parameter vector")
  }
  if(length(extra_estimated_params)>0){
    assert_that(all(params_data$name[(2*n_env_vars)+c(1:n_extra_params)] == extra_estimated_params),
                msg = "Initial parameter vector must end with additional estimated parameters")
  }

  return(NULL)
}

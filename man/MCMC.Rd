% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{MCMC}
\alias{MCMC}
\title{MCMC}
\usage{
MCMC(
  log_params_ini = c(),
  input_data = list(),
  obs_sero_data = NULL,
  obs_case_data = NULL,
  filename_prefix = "Chain",
  Niter = 1,
  type = NULL,
  mode_start = 0,
  prior_settings = list(type = "zero"),
  dt = 1,
  n_reps = 1,
  enviro_data = NULL,
  R0_fixed_values = NULL,
  p_severe_inf = 0.12,
  p_death_severe_inf = 0.39,
  add_values = list(vaccine_efficacy = 1, p_rep_severe = 1, p_rep_death = 1, m_FOI_Brazil
    = 1),
  deterministic = FALSE,
  mode_parallel = "none",
  cluster = NULL
)
}
\arguments{
\item{log_params_ini}{Initial values of parameters to be varied. These should always be the log() values of the
actual parameters, ordered as follows:
1) Parameters controlling the value of spillover force of infection FOI, either a) a number of FOI values equal
to the total number of regions to be considered or b) a number of environmental coefficients used to calculate
FOI values from environmental covariates equal to the number of environmental covariates listed in the
enviro_data frame. Values should be in alphabetical order by region in case (a) or in the order of the columns
in the environmental data frame in case (b).
2) If the basic reproduction number for human-human transmission R0 is to be varied (i.e. type is set to
"FOI+R0" or "FOI+R0 enviro"), parameters controlling the value of R0, either a) a number of R0 values equal to
the total number of regions to be considered or b) a number of environmental coefficients used to calculate R0
values from environmental covariates equal to the number of environmental covariates listed in the enviro_data
frame. Values should be in alphabetical order by region in case (a) or in the order of the columns in the
environmental data frame in case (b).
3) Values of the additional parameters (vaccine efficacy vaccine_efficacy, severe case reporting probability
p_rep_severe and fatal case reporting probability p_rep_death) if these are to be varied, in the order
vaccine_efficacy->p_rep_severe->p_rep_death. If these parameters are to be varied, the values separately
supplied to this function (see below) should be set to NULL, the default.}

\item{input_data}{List of population and vaccination data for multiple regions (created using data input creation
code and usually loaded from RDS file)}

\item{obs_sero_data}{Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
positives}

\item{obs_case_data}{Annual reported case/death data for comparison, by region and year, in format no.
cases/no. deaths}

\item{filename_prefix}{Prefix of names for output files}

\item{Niter}{Total number of steps to run}

\item{type}{Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination
If mode_start=0, only vaccinated individuals
If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
If mode_start=3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)}

\item{prior_settings}{List containing settings for priors: must contain text named "type":
If type = "zero", prior probability is always zero
If type = "flat", prior probability is zero if log parameter values in designated ranges log_params_min and log_params_max,
 -Inf otherwise; log_params_min and log_params_max included in prior_settings as vectors of same length as log_params_ini
If type = "exp", prior probability is given by dexp calculation on FOI/R0 values
If type = "norm", prior probability is given by dnorm calculation on parameter values with settings based on vectors of values
 in prior_settings; norm_params_mean and norm_params_sd (mean and standard deviation values applied to log FOI/R0
 parameters and to actual values of additional parameters) + R0_mean + R0_sd (mean + standard deviation of computed R0, single values)}

\item{dt}{time increment in days (must be 1 or 5)}

\item{n_reps}{Number of times to repeat calculations to get average likelihood at each step}

\item{enviro_data}{Data frame containing values of environmental covariates; set to NULL if not in use}

\item{R0_fixed_values}{Values of R0 to use if only FOI is subject to fitting (i.e. type set to "FOI" or "FOI
enviro"); set to NULL if not in use}

\item{p_severe_inf}{Probability of an infection being severe}

\item{p_death_severe_inf}{Probability of a severe infection resulting in death}

\item{add_values}{List containing values of additional inputs or NULL to indicate they are being varied as parameters:
vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (set to NULL if being varied as a parameter)}

\item{deterministic}{TRUE/FALSE - set model to run in deterministic mode if TRUE}

\item{mode_parallel}{Set mode for parallelization, if any:
If mode_parallel="none", no parallelization
If mode_parallel="pars_multi", all regions run in parallel for same time period with same output type
If mode_parallel="clusterMap", all regions run in parallel with different time periods and output types}

\item{cluster}{Cluster of threads to use if mode_parallel="clusterMap"
'}
}
\description{
Combined MCMC Multi-Region - series of MCMC steps for one or more regions
}
\details{
This is the master function for running a Markov chain to optimize the parameters of the yellow fever
model based on the calculated likelihood of observing supplied data given a particular set of parameters.
}

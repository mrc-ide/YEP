% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{MCMC}
\alias{MCMC}
\title{MCMC}
\usage{
MCMC(
  params_data = data.frame(name = "FOI_var1", initial = 1, max = Inf, min = -Inf, mean =
    0, sd = 1, estimate = TRUE),
  input_data = list(),
  obs_sero_data = NULL,
  obs_case_data = NULL,
  filename_prefix = "Chain",
  Niter = 1,
  mode_start = 1,
  time_inc = 1,
  n_reps = 1,
  enviro_data_const = list(),
  enviro_data_var = list(),
  deterministic = FALSE,
  mode_time = 1,
  mode_parallel = FALSE,
  cluster = NULL
)
}
\arguments{
\item{params_data}{#Data frame of parameter information containing names, initial values, maximum and minimum values,
mean and standard deviation (for prior calculation) and flag indicating whether parameter estimated or fixed\cr
Parameters to include: coefficients of environmental covariates to calculate FOI_spillover and R0, reported
vaccination effectiveness, probability of severe case reporting, probability of fatal case reporting, Brazil
FOI_spillover multiplier, FOI_spillover and R0 (latter two never estimated, included only for priors)\cr
TBA - instructions}

\item{input_data}{List of population and vaccination data for multiple regions (created using data input creation
code and usually loaded from RDS file)}

\item{obs_sero_data}{Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
positives (TBA - instructions)}

\item{obs_case_data}{Annual reported case/death data for comparison, by region and year, in format no.
cases/no. deaths (TBA - instructions)}

\item{filename_prefix}{Prefix of output RDS file name, e.g. "Chain.Rds"}

\item{Niter}{Total number of iterations to run}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination \cr
If mode_start = 0, only vaccinated individuals \cr
If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr}

\item{time_inc}{time increment in days (must be 1 or 5)}

\item{n_reps}{Number of times to repeat calculations to get average likelihood at each iteration}

\item{enviro_data_const}{Data frame of values of constant environmental covariates (columns) by region (rows)}

\item{enviro_data_var}{List containing time-varying environmental covariate data:\cr
regions: Vector of region labels\cr
env_vars: Vector of covariate names\cr
values: Array of covariate values with dimensions (number of covariates, number of regions, number of time points).
Number of time points must be correct for mode_time setting.\cr}

\item{deterministic}{TRUE/FALSE - set model to run in deterministic mode if TRUE}

\item{mode_time}{Type of time dependence of FOI_spillover and R0 to be used: \cr
If mode_time = 0, no time variation (constant values)\cr
If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)}

\item{mode_parallel}{TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed}

\item{cluster}{Cluster of threads to use if mode_parallel = TRUE
'}
}
\description{
Combined MCMC Multi-Region - series of MCMC iterations for one or more regions
}
\details{
This is the master function for running a Markov chain to optimize the parameters of the yellow fever
model based on the calculated likelihood of observing supplied data given a particular set of parameters.
}

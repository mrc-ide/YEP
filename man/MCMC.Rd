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
  log_params_min = c(),
  log_params_max = c(),
  mode_start = 0,
  prior_type = "zero",
  dt = 1,
  enviro_data = NULL,
  R0_fixed_values = NULL,
  vaccine_efficacy = NULL,
  p_rep_severe = NULL,
  p_rep_death = NULL,
  m_FOI_Brazil = 1,
  cluster = NULL
)
}
\arguments{
\item{log_params_ini}{Initial values of parameters to be varied These should always be the log() values of the
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

\item{log_params_min}{Lower limits of varied parameter values if specified}

\item{log_params_max}{Upper limits of varied parameter values if specified}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination
If mode_start = 0, only vaccinated individuals
If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity}

\item{prior_type}{Text indicating which type of calculation to use for prior probability
If prior_type = "zero", prior probability is always zero
If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
If prior_type = "norm", prior probability is given by dnorm calculation on parameter values}

\item{dt}{time increment in days (must be 1 or 5)}

\item{enviro_data}{Data frame containing values of environmental covariates; set to NULL if not in use}

\item{R0_fixed_values}{Values of R0 to use if only FOI is subject to fitting (i.e. type set to "FOI" or "FOI
enviro"); set to NULL if not in use}

\item{vaccine_efficacy}{Vaccine efficacy (set to NULL if being varied as a parameter)}

\item{p_rep_severe}{Probability of observation of severe infection (set to NULL if being varied as a parameter)}

\item{p_rep_death}{Probability of observation of death (set to NULL if being varied as a parameter)}

\item{m_FOI_Brazil}{Multiplier of spillover FOI for Brazil regions (set to NULL if being varied as a parameter)}

\item{cluster}{Cluster of threads to use for multithreading; set to NULL if not using multithreading
'}
}
\description{
Combined MCMC Multi-Region - series of MCMC steps for one or more regions
}
\details{
This is the master function for running a Markov chain to optimize the parameters of the yellow fever
model based on the calculated likelihood of observing supplied data given a particular set of parameters.
}

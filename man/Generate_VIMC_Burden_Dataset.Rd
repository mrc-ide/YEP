% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputs_generate.R
\name{Generate_VIMC_Burden_Dataset}
\alias{Generate_VIMC_Burden_Dataset}
\title{Generate_VIMC_Burden_Dataset}
\usage{
Generate_VIMC_Burden_Dataset(
  input_data = list(),
  FOI_values = c(),
  R0_values = c(),
  template = NULL,
  vaccine_efficacy = 1,
  p_severe_inf = 0.12,
  p_death_severe_inf = 0.39,
  YLD_per_case = 0.006486,
  mode_start = 1,
  start_SEIRV = NULL,
  dt = 1,
  n_reps = 1,
  deterministic = FALSE,
  mode_parallel = FALSE,
  cluster = NULL,
  seed = NULL
)
}
\arguments{
\item{input_data}{List of population and vaccination data for multiple regions}

\item{FOI_values}{Values for each region of the force of infection due to spillover from sylvatic reservoir}

\item{R0_values}{Values for each region of the basic reproduction number for human-human transmission}

\item{template}{Burden data in VIMC format, with regions, years, minimum and maximum age, and life expectancy for each line}

\item{vaccine_efficacy}{Fractional vaccine efficacy}

\item{p_severe_inf}{Probability of an infection being severe}

\item{p_death_severe_inf}{Probability of a severe infection resulting in death}

\item{YLD_per_case}{TBA}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination
If mode_start=0, only vaccinated individuals
If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
If mode_start=2, use SEIRV input in list from previous run(s)
If mode_start=3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)}

\item{start_SEIRV}{SEIRV data from end of a previous run to use as input (list of datasets, one per region)}

\item{dt}{Time increment in days to use in model (should be either 1.0 or 5.0 days)}

\item{n_reps}{number of stochastic repetitions}

\item{deterministic}{TRUE/FALSE - set model to run in deterministic mode if TRUE}

\item{mode_parallel}{TRUE/FALSE - set model to run in parallel using cluster if TRUE}

\item{cluster}{Cluster of threads to use if mode_parallel=TRUE}

\item{seed}{Optional random seed value to set before running each region for stochastic normalization; set to NULL
  to omit; will not work if mode_parallel is not set to "none".
'}
}
\description{
Generate annual burden data in format used by VIMC
}
\details{
This function is used to generate annual burden data in the format used by the Vaccine Impact
  Modelling Consortium (VIMC) [TBA]
}

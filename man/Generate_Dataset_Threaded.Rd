% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel.R
\name{Generate_Dataset_Threaded}
\alias{Generate_Dataset_Threaded}
\title{Generate_Dataset_Threaded}
\usage{
Generate_Dataset_Threaded(
  input_data,
  FOI_values,
  R0_values,
  obs_sero_data,
  obs_case_data,
  vaccine_efficacy,
  p_rep_severe,
  p_rep_death,
  mode_start,
  dt,
  cluster
)
}
\arguments{
\item{input_data}{List of population and vaccination data for multiple regions}

\item{FOI_values}{Values for each region of the force of infection due to spillover from sylvatic reservoir}

\item{R0_values}{Values for each region of the basic reproduction number for human-human transmission}

\item{obs_sero_data}{Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
positives}

\item{obs_case_data}{Annual reported case/death data for comparison, by region and year, in format no. cases/no.
deaths}

\item{vaccine_efficacy}{Fractional vaccine efficacy}

\item{p_rep_severe}{Probability of reporting of a severe but non-fatal infection}

\item{p_rep_death}{Probability of reporting of a fatal infection}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination
If mode_start=0, only vaccinated individuals
If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
If mode_start=2, use SEIRV input in list from previous run(s)}

\item{dt}{Time increment in days to use in model (should be either 1.0 or 5.0 days)}

\item{cluster}{Cluster of threads to use for multithread run
'}
}
\description{
Generate serological and/or annual case/death data using multithreading
}
\details{
Version of Generate_Dataset function for multithreading
}

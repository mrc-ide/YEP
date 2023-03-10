# R file for functions used for multithread dataset generation in YEP package
#-------------------------------------------------------------------------------
#' @title Model_Run_Thread
#'
#' @description Run SEIRV model for one region when dataset being generated using multithreading
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs specific data for dataset being
#' generated in parallel.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years for which to output SEIRV data
#' @param flag_case TBA
#' @param flag_sero TBA
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRVC input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
Model_Run_Thread <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),years_data=c(),flag_case=0,flag_sero=0,
                             year0=1940,mode_start=0,vaccine_efficacy=1.0,dt=1.0){

  model_output=Model_Run(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,vaccine_efficacy,dt=dt)

  model_data=list(years=years_data)
  if(flag_case==1){
    model_data$infs=rep(NA,length(years_data))
    i=0
    for(year in years_data){
      i=i+1
      model_data$infs[i]=sum(model_output$C[model_output$year==year,])
    }
  }
  if(flag_sero==1){
    N_age=dim(vacc_data)[2]
    blank=array(NA,dim=c(length(years_data),N_age))
    model_data$S=model_data$E=model_data$I=model_data$R=model_data$V=blank
    i=0
    for(year in years_data){
      i=i+1
      for(j in 1:N_age){
        model_data$S[i,j]=sum(model_output$S[model_output$year==year,j])
        model_data$E[i,j]=sum(model_output$E[model_output$year==year,j])
        model_data$I[i,j]=sum(model_output$I[model_output$year==year,j])
        model_data$R[i,j]=sum(model_output$R[model_output$year==year,j])
        model_data$V[i,j]=sum(model_output$V[model_output$year==year,j])
      }
    }
  }
  model_output=NULL

  return(model_data)
}
#-------------------------------------------------------------------------------
#' @title Generate_Dataset_Threaded
#'
#' @description Generate serological and/or annual case/death data using multithreading
#'
#' @details Version of Generate_Dataset function for multithreading
#'
#' @param input_data List of population and vaccination data for multiple regions
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param cluster Cluster of threads to use for multithread run
#' '
#' @export
#'
Generate_Dataset_Threaded <- function(input_data,FOI_values,R0_values,obs_sero_data,obs_case_data,
                                      vaccine_efficacy,p_rep_severe,p_rep_death,mode_start,dt,cluster){

  assert_that(is.null(cluster)==FALSE)
  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                            " (see https://mrc-ide.github.io/YellowFeverDynamics/articles/CGuideAInputs.html )"))
  assert_that(any(is.null(obs_sero_data)==FALSE,is.null(obs_case_data)==FALSE),
              msg="Need at least one of obs_sero_data or obs_case_data")
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0,msg="Vaccine efficacy must be between 0 and 1")
  if(is.null(obs_case_data)==FALSE){
    assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0,
                msg="Severe infection reporting probability must be between 0 and 1")
    assert_that(p_rep_death >=0.0 && p_rep_death <=1.0,
                msg="Fatal infection reporting probability must be between 0 and 1")
  }

  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions")

  if(is.null(input_data$flag_sero)){
    input_data=input_data_process(input_data,obs_sero_data,obs_case_data)
  }

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(obs_sero_data)==FALSE){
    blank=rep(0,nrow(obs_sero_data))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  } else {model_sero_data=NULL}

  if(is.null(obs_case_data)==FALSE){
    model_case_values=model_death_values=rep(0,nrow(obs_case_data))
  } else {model_case_values=model_death_values=NA}


  #Model all regions and save relevant output data
  vacc_data_subsets=pop_data_subsets=years_data_sets=list()
  for(n_region in 1:n_regions){
    vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
    pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
    years_data_sets[[n_region]]=c(input_data$year_data_begin[n_region]:input_data$year_end[n_region])
  }

  model_data <- clusterMap(cluster,Model_Run_Thread,FOI=FOI_values,R0=R0_values,vacc_data=vacc_data_subsets,
                           pop_data=pop_data_subsets,years_data=years_data_sets,flag_case=input_data$flag_case,
                           flag_sero=input_data$flag_sero,
                           MoreArgs=list(year0=input_data$years_labels[1],mode_start=mode_start,
                                         vaccine_efficacy=vaccine_efficacy,dt=dt))

  for(n_region in 1:n_regions){

    #Get model output data for region
    model_data_subset=model_data[[n_region]]

    #Compile case data if needed
    if(input_data$flag_case[n_region]==1){
      case_line_list=input_data$case_line_list[[n_region]]
      years_outbreak=obs_case_data$year[case_line_list]
      n_years_outbreak=length(case_line_list)

      rep_cases=rep_deaths=rep(0,n_years_outbreak)
      for(n_year in 1:n_years_outbreak){
        year=years_outbreak[n_year]
        infs=model_data_subset$infs[model_data_subset$years==year]
        severe_infs=rbinom(1,floor(infs),p_severe_inf)
        deaths=rbinom(1,severe_infs,p_death_severe_inf)
        rep_deaths[n_year]=rbinom(1,deaths,p_rep_death)
        rep_cases[n_year]=rep_deaths[n_year]+rbinom(1,severe_infs-deaths,p_rep_severe)
      }

      model_case_values[case_line_list]=model_case_values[case_line_list]+rep_cases
      model_death_values[case_line_list]=model_death_values[case_line_list]+rep_deaths
    }

    #Compile seroprevalence data if necessary
    if(input_data$flag_sero[n_region]==1){
      sero_line_list=input_data$sero_line_list[[n_region]]
      sero_results=sero_calculate2(obs_sero_data[sero_line_list,],model_data_subset)
      model_sero_data$samples[sero_line_list]=model_sero_data$samples[sero_line_list]+sero_results$samples
      model_sero_data$positives[sero_line_list]=model_sero_data$positives[sero_line_list]+sero_results$positives
    }
  }
  model_output<-NULL

  if(any(input_data$flag_sero==1)){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}

  return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
              model_death_values=model_death_values))
}

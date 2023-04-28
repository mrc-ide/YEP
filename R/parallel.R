# R file for functions used for multithread dataset generation in YEP package
#-------------------------------------------------------------------------------
#' @title Model_Run_Threaded
#'
#' @description Run SEIRV model for one region when dataset being generated using multithreading
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs specific data for dataset
#' being generated in parallel.
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
Model_Run_Threaded <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),years_data=c(),flag_case=0,
                             flag_sero=0,year0=1940,mode_start=0,vaccine_efficacy=1.0,dt=1.0){

  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,vaccine_efficacy,dt=dt)

  x <- SEIRV_Model$new(FOI_spillover=pars$FOI_spillover,R0=pars$R0,vacc_rate_annual=pars$vacc_rate_annual,
                       steps_inc=pars$steps_inc,steps_lat=pars$steps_lat,steps_inf=pars$steps_inf,
                       Cas0=pars$Cas0,Exp0=pars$Exp0,Inf0=pars$Inf0,N_age=pars$N_age,Rec0=pars$Rec0,Sus0=pars$Sus0,
                       Vac0=pars$Vac0,dP1_all=pars$dP1_all,dP2_all=pars$dP2_all,n_years=pars$n_years,
                       year0=pars$year0,vaccine_efficacy=pars$vaccine_efficacy,dt=pars$dt)

  n_nv=4 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  t_pts_all=c(1:((max(years_data)+1-year0)*(365/dt))) #All output time points
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  n_steps=length(t_pts_all) #Total number of output time points
  step0=(years_data[1]-year0)*(365/dt) #Step at which data starts being saved for final output
  t_pts_out=n_steps-step0 #Number of time points in final output data

  x_res <- x$run(n_steps)
  t_pts=c((step0+1):n_steps)
  if(step0==0){x_res[1,3]=year0}

  model_output=list(year=x_res[t_pts,3])
  if(flag_sero==1){
    model_output$S=array(x_res[t_pts,c((1+n_nv):(N_age+n_nv))],dim=c(t_pts_out,N_age))
    model_output$E=array(x_res[t_pts,c((N_age+1+n_nv):((2*N_age)+n_nv))],dim=c(t_pts_out,N_age))
    model_output$I=array(x_res[t_pts,c(((2*N_age)+1+n_nv):((3*N_age)+n_nv))],dim=c(t_pts_out,N_age))
    model_output$R=array(x_res[t_pts,c(((3*N_age)+1+n_nv):((4*N_age)+n_nv))],dim=c(t_pts_out,N_age))
    model_output$V=array(x_res[t_pts,c(((4*N_age)+1+n_nv):((5*N_age)+n_nv))],dim=c(t_pts_out,N_age))
  }
  if(flag_case==1) {
    model_output$C=array(x_res[t_pts,c(((5*N_age)+1+n_nv):((6*N_age)+n_nv))],dim=c(t_pts_out,N_age))
  }

  model_data=list(years=years_data)
  t_pts=length(model_output$year)
  n_years=length(years_data)
  if(flag_case==1){
    model_data$infs=rep(NA,n_years)
    for(n_year in c(1:n_years)){
      pts=c(1:t_pts)[model_output$year==years_data[n_year]]
      model_data$infs[n_year]=sum(model_output$C[pts,])
    }
  }
  if(flag_sero==1){
    N_age=dim(vacc_data)[2]
    blank=array(NA,dim=c(n_years,N_age))
    model_data$S=model_data$E=model_data$I=model_data$R=model_data$V=blank
    for(n_year in c(1:n_years)){
      year=years_data[n_year]
      pts=c(1:t_pts)[model_output$year==year]
      for(j in 1:N_age){
        model_data$S[n_year,j]=sum(model_output$S[pts,j])
        model_data$E[n_year,j]=sum(model_output$E[pts,j])
        model_data$I[n_year,j]=sum(model_output$I[pts,j])
        model_data$R[n_year,j]=sum(model_output$R[pts,j])
        model_data$V[n_year,j]=sum(model_output$V[pts,j])
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
  assert_that(input_data_check(input_data),msg="Input data must be in standard format")
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

  model_data <- clusterMap(cl=cluster,fun=Model_Run_Threaded,FOI_spillover=FOI_values,R0=R0_values,vacc_data=vacc_data_subsets,
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
        infs=sum(model_data_subset$infs[model_data_subset$years==years_outbreak[n_year]])
        severe_infs=infs*p_severe_inf
        deaths=severe_infs*p_death_severe_inf
        rep_deaths[n_year]=round(deaths*p_rep_death)
        rep_cases[n_year]=rep_deaths[n_year]+round((severe_infs-deaths)*p_rep_severe)
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

  if(any(input_data$flag_sero>0)){
    model_sero_data$sero=model_sero_data$positives/model_sero_data$samples
  }

  #model_data=NULL

  return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
              model_death_values=model_death_values))
}

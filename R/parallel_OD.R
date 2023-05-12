#-------------------------------------------------------------------------------
#' @title Model_Run_Multi_Region
#'
#' @description Run SEIRV model for multiple regions
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for multiple regions over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values.
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#' @param R0 Vector of values of basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each region and age group by year
#' @param pop_data Population in each region and age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param output_type Type of data to output: "full" = SEIRVC, "case" = C only, "sero" = SEIRV only
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Multi_Region <- function(FOI_spillover = c(),R0 = c(),vacc_data = list(),pop_data = list(),
                                   years_data = c(1940:1941),start_SEIRV = list(), output_type = "full", year0 = 1940,
                                   mode_start = 0,vaccine_efficacy = 1.0,dt = 1.0,n_particles = 1, n_threads = 1,
                                   deterministic = FALSE) {

  #TODO Add assert_that functions

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,1,]) #Number of age groups
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data
  dimensions=c(N_age,n_particles,t_pts_out)

  n_regions=length(FOI_spillover)
  pars=list()
  for(i in 1:n_regions){
    pars[[i]]=parameter_setup(FOI_spillover[i],R0[i],vacc_data[i,,],pop_data[i,,],year0,years_data,mode_start,
                              vaccine_efficacy,start_SEIRV[[i]],dt)
  }

  x <- SEIRV_Model$new(pars,time = 1, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic,
                       pars_multi = TRUE)

  x_res <- array(NA, dim = c(n_data_pts, n_particles*n_regions, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles*n_regions)}

  output=list()
  for(i in 1:n_regions){
    n_p2=c(1:n_particles)+(n_particles*(i-1))
    if(output_type=="full"){
      dimensions=c(N_age,n_particles,t_pts_out)
      output[[i]]=list(day=x_res[1,n_p2[1],],year=x_res[2,n_p2[1],])
      output[[i]]$FOI_total=array(x_res[3,n_p2,]/dt,dim=c(n_particles,t_pts_out))
      output[[i]]$S=array(x_res[c((1+n_nv):(N_age+n_nv)),n_p2,],dim=dimensions)
      output[[i]]$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p2,],dim=dimensions)
      output[[i]]$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),n_p2,],dim=dimensions)
      output[[i]]$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),n_p2,],dim=dimensions)
      output[[i]]$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),n_p2,],dim=dimensions)
      output[[i]]$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p2,],dim=dimensions)
    } else {
      n_years=length(years_data)
      dimensions=c(N_age,n_particles,n_years)
      output[[i]]=list(year=years_data)
      if(output_type=="case+sero" || output_type=="sero"){
        output[[i]]$V=output[[i]]$R=output[[i]]$I=output[[i]]$E=output[[i]]$S=array(0,dim=dimensions)
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,n_p2[1],]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output[[i]]$S[,n_p,n_year]=rowMeans(x_res[c((1+n_nv):(N_age+n_nv)),n_p2[n_p],pts])
            output[[i]]$E[,n_p,n_year]=rowMeans(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p2[n_p],pts])
            output[[i]]$I[,n_p,n_year]=rowMeans(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),n_p2[n_p],pts])
            output[[i]]$R[,n_p,n_year]=rowMeans(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),n_p2[n_p],pts])
            output[[i]]$V[,n_p,n_year]=rowMeans(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),n_p2[n_p],pts])
          }
        }
      }
      if(output_type=="case+sero" || output_type=="case"){
        output[[i]]$C=array(0,dim=c(n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,n_p2[1],]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output[[i]]$C[n_p,n_year]=sum(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p2[n_p],pts])
          }
        }
      }
    }
  }

  return(output)
}

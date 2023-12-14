#' #' @title Model_Run2
#' #'
#' #' @description Run SEIRV model for single region (alternative version with amendments to increase speed in some cases)
#' #'
#' #' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' #' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' #' values, infection numbers and/or total force of infection values.
#' #'
#' #' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' #' @param R0 Basic reproduction number for urban spread of infection
#' #' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy=1) by age group and year
#' #' @param pop_data Population by age group and year
#' #' @param years_data Incremental vector of years denoting years for which to save data
#' #' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' #' @param output_type Type of data to output:
#' #'   "full" = SEIRVC + FOI for all steps and ages
#' #'   "case" = annual total new infections (C) summed across all ages
#' #'   "sero" = annual SEIRV
#' #'   "case+sero" = annual SEIRVC, cases summed across all ages
#' #'   "case_alt" = annual total new infections not combined by age
#' #'   "case_alt2" = total new infections combined by age for all steps
#' #' @param year0 First year in population/vaccination data
#' #' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#' #'  If mode_start=0, only vaccinated individuals
#' #'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
#' #'  If mode_start=2, use SEIRV input in list from previous run(s)
#' #'  If mode_start=3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' #' @param vaccine_efficacy Proportional vaccine efficacy
#' #' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' #' @param n_particles number of particles to use
#' #' @param n_threads number of threads to use
#' #' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' #' '
#' #' @export
#' #'
#' Model_Run2 <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
#'                       start_SEIRV = list(), output_type = "full", year0 = 1940, mode_start = 0,
#'                       vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {
#' 
#'   #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
#'   assert_that(n_particles<=20,msg="Number of particles must be 20 or less")
#' 
#'   n_nv=3 #Number of non-vector outputs
#'   N_age=length(pop_data[1,]) #Number of age groups
#'   #n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
#'   step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
#'   step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
#'   t_pts_out=step_end-step_begin+1 #Number of time points in final output data
#' 
#'   x <- SEIRV_Model$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
#'                                             vaccine_efficacy,start_SEIRV,dt),
#'                        time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)
#' 
#'   if(output_type %in% c("full","sero","case+sero")){
#'     x_res=array(NA,dim=c(n_nv+(6*N_age),n_particles,t_pts_out))
#'     for(step in step_begin:step_end){
#'       x_res[,,step-step_begin+1] <- x$run(step)
#'     }
#'   }else{
#'     indices_out=c(1,2,x$info()$index$C)
#'     x_res=array(NA,dim=c(length(indices_out),n_particles,t_pts_out))
#'     for(step in step_begin:step_end){
#'       x_res[,,step-step_begin+1] <- x$run(step)[indices_out,]
#'     }
#'   }
#'   if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}
#' 
#'   if(output_type=="full"){
#'     dimensions=c(N_age,n_particles,t_pts_out)
#'     output_data=list(day=x_res[1,1,],year=x_res[2,1,])
#'     output_data$FOI_total=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
#'     output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
#'     output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
#'     output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dimensions)
#'     output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dimensions)
#'     output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dimensions)
#'     output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dimensions)
#'   } else {
#'     if(output_type=="case_alt2"){
#'       output_data=list(day=x_res[1,1,],year=x_res[2,1,])
#'       output_data$C=array(0,dim=c(n_particles,t_pts_out))
#'       pts2=c(3:(N_age+2))
#'       for(pt in 1:t_pts_out){
#'         for(n_p in 1:n_particles){
#'           output_data$C[n_p,pt]=sum(x_res[pts2,n_p,pt])
#'         }
#'       }
#'     }  else {
#'       n_years=length(years_data)
#'       output_data=list(year=years_data)
#'       if(output_type=="case+sero" || output_type=="sero"){
#'         output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=c(N_age,n_particles,n_years))
#'         for(n_year in 1:n_years){
#'           pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
#'           for(n_p in 1:n_particles){
#'             output_data$S[,n_p,n_year]=rowMeans(x_res[c((1+n_nv):(N_age+n_nv)),n_p,pts])
#'             output_data$E[,n_p,n_year]=rowMeans(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p,pts])
#'             output_data$I[,n_p,n_year]=rowMeans(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),n_p,pts])
#'             output_data$R[,n_p,n_year]=rowMeans(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),n_p,pts])
#'             output_data$V[,n_p,n_year]=rowMeans(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),n_p,pts])
#'           }
#'         }
#'       }
#'       if(output_type=="case+sero" || output_type=="case"){
#'         output_data$C=array(0,dim=c(n_particles,n_years))
#'         if(output_type=="case"){pts2=c(3:(N_age+2))}else{pts2=c(((5*N_age)+1+n_nv):((6*N_age)+n_nv))}
#'         for(n_year in 1:n_years){
#'           pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
#'           for(n_p in 1:n_particles){
#'             output_data$C[n_p,n_year]=sum(x_res[pts2,n_p,pts])
#'           }
#'         }
#'       }
#'       if(output_type=="case_alt"){
#'         output_data$C=array(0,dim=c(N_age,n_particles,n_years))
#'         pts2=c(3:(N_age+2))
#'         for(n_year in 1:n_years){
#'           pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
#'           for(n_p in 1:n_particles){
#'             output_data$C[,n_p,n_year]=rowSums(x_res[pts2,n_p,pts])
#'           }
#'         }
#'       }
#'     }
#'   }
#'   return(output_data)
#' }
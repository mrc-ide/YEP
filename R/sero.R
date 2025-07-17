# R file for functions relating to serological data in YEP package
#-------------------------------------------------------------------------------
#TODO - Update for new model versions
#' @title sero_calculate
#'
#' @description Calculate seroprevalence in unvaccinated people from modelled data for one or more years and one age
#' range
#'
#' @details Takes in information on minimum and maximum ages of desired range, year(s) for which to calculate
#' seroprevalence, factor representing proportion of patients with unknown vaccine status, and SEIRV model output
#' data, and calculates seroprevalence in unvaccinated people in specified age range for specified year(s).
#'
#' @param age_min = Minimum age of age group
#' @param age_max = Maximum age of age group
#' @param years = Years for which to calculate average annual seroprevalence
#' @param vc_factor = Proportion of patients tested for whom vaccine status unknown
#' @param data = Output of Model_Run
#' @param n_p = Particle to select from data
#' '
#' @export
#'
sero_calculate <- function(age_min=0,age_max=101,years=NULL,vc_factor=0,data=list(),n_p=1){

  assert_that(age_min>=0,msg="Minimum age must be equal to or greater than 0")
  assert_that(age_max>age_min,msg="Maximum age must be greater than minimum age")
  assert_that(is.null(years)==FALSE,msg="Years in which to calculate seroprevalence must be specified")
  assert_that(between(vc_factor,0.0,1.0),msg = "vc_factor must be between 0-1")
  assert_that(is.null(data$S)==FALSE) #TODO - Improve check on SEIRV data
  assert_that(n_p %in% c(1:dim(data$S)[2]),
              msg="Selected particle number must be a positive integer equal to or less than number of particles")

  ages=c((age_min+1):age_max)
  sero_values=rep(0,length(years))

  for(i in 1:length(years)){
    n_t=which(data$year %in% years[i])
    if(vc_factor==0){
      samples=data$S[ages,n_p,n_t]+data$E[ages,n_p,n_t]+data$I[ages,n_p,n_t]+data$R[ages,n_p,n_t]
      positives=data$R[ages,n_p,n_t]
      sero_values[i]=sum(positives)/sum(samples)
    } else {
      if(vc_factor==1){
        samples=data$S[ages,n_p,n_t]+data$E[ages,n_p,n_t]+data$I[ages,n_p,n_t]+data$R[ages,n_p,n_t]+data$V[ages,n_p,n_t]
        positives=data$R[ages,n_p,n_t]+data$V[ages,n_p,n_t]
        sero_values[i]=sum(positives)/sum(samples)
      } else {
        samples=data$S[ages,n_p,n_t]+data$E[ages,n_p,n_t]+data$I[ages,n_p,n_t]+data$R[ages,n_p,n_t]
        positives=data$R[ages,n_p,n_t]
        sero_values[i]=((1.0-vc_factor)*sum(positives))/sum(samples)
        samples=samples+data$V[ages,n_p,n_t]
        positives=positives+data$V[ages,n_p,n_t]
        sero_values[i]=sero_values[i]+((vc_factor*sum(positives))/sum(samples))
      }
    }
  }

  return(sero_values)
}
#-------------------------------------------------------------------------------
#' @title sero_calculate2_alt
#'
#' @description Calculate number of "samples" and number of "positives" from modelled data for specified age range(s)
#' and year(s)
#'
#' @details Takes in information on minimum and maximum ages of desired range(s), year(s) for which to calculate
#' number of "samples" (people eligible for testing) and "positives" (people who would test positive), plus vc_factor
#' (proportion of people for whom vaccination status unknown)
#'
#' @param sero_data Data frame containing years, minimum and maximum ages, and values of vc_factor (proportion of
#' people for whom vaccination status unknown)
#' @param model_data Annual cumulative SEIRV output of Model_Run2
#' @param n_region Region number to select from model_data
#' @param n_p Particle to select from model_data
#' '
#' @export
#'
sero_calculate2_alt <- function(sero_data=list(),model_data=list(),n_region=1,n_p=1){
  assert_that(is.data.frame(sero_data))
  assert_that(is.list(model_data))
  assert_that(n_p<=dim(model_data$SEIR_annual)[3],msg="Specified particle number is unavailable")

  nrows=nrow(sero_data)
  output_frame=data.frame(samples=rep(NA,nrows),positives=rep(NA,nrows))

  for(i in 1:nrows){
    ages=c((sero_data$age_min[i]+1):sero_data$age_max[i])
    year=sero_data$year[i]
    vc_factor=sero_data$vc_factor[i]
    n_t=which(model_data$year==year)
    SEIR_sum=sum(model_data$SEIR_annual[n_region,ages,n_p,n_t])
    R_sum=sum(model_data$R_annual[n_region,ages,n_p,n_t])
    if(vc_factor>0){
      V_sum=sum(model_data$V_annual[n_region,ages,n_p,n_t])
      if(vc_factor==1){
        samples=SEIR_sum+V_sum
        positives=R_sum+V_sum
      } else {
        samples=SEIR_sum
        T_sum=SEIR_sum+V_sum
        positives=((1.0-vc_factor)*R_sum)+(vc_factor*(SEIR_sum/T_sum)*(R_sum+V_sum))
      }
    } else {
      samples=SEIR_sum
      positives=R_sum
    }
    output_frame$samples[i]=samples
    output_frame$positives[i]=positives
  }

  return(output_frame)
}
#-------------------------------------------------------------------------------
#' @title sero_data_compare
#'
#' @description Take seroprevalence results from dataset generation functions, compare comparison with
#' observed/template seroprevalence data and calculate likelihood
#'
#' @details [TBA]
#'
#' @param model_sero_values Seroprevalence values from dataset generation function (no. positives/no. samples)
#' @param obs_sero_data Seroprevalence data for comparison, by year and age group, in format
#'   no. samples/no. positives
#' '
#' @export
#'
sero_data_compare <- function(model_sero_values=c(),obs_sero_data=list()){

  #TODO - Add assert_that functions

  #Eliminate values which can cause -Inf results
  model_sero_values[model_sero_values>1.0]=1.0
  model_sero_values[is.infinite(model_sero_values)]=0.0
  sero_rev=1.0-model_sero_values
  sero_rev[sero_rev<0.0]=0.0

  sero_like_values=lgamma(obs_sero_data$samples+1)-lgamma(obs_sero_data$positives+1)-
    lgamma(obs_sero_data$samples-obs_sero_data$positives+1)+obs_sero_data$positives*log(model_sero_values)+
    (obs_sero_data$samples-obs_sero_data$positives)*log(sero_rev)

  return(sero_like_values)
}

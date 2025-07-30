# File for functions used to estimate parameters via Monty
# TODO:
# -Add check function(s) and/or assert_that checks in existing functions
# -Make separate function for new region indexing?
#-------------------------------------------------------------------------------
#' @title pars_var_setup
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param n_env_vars TBA
#' @param vars_extra_names TBA
#'
#' @export
#'
pars_var_setup <- function(n_env_vars=5,vars_extra_names=c("p_rep_severe","p_rep_death")){

  assertthat::assert_that(all(vars_extra_names %in% extra_param_names))
  n_extra = length(vars_extra_names)
  n_rows = n_extra + (2*n_env_vars)
  pars_var=data.frame(name=c(vars_extra_names,paste0("log_FOI_coeffs[",c(1:n_env_vars),"]"),
                             paste0("log_R0_coeffs[",c(1:n_env_vars),"]")),
                      initial=c(rep(1,n_extra),rep(-15,n_env_vars),rep(-3,n_env_vars)),
                      max=c(rep(1,n_extra),rep(-10,n_env_vars),rep(1,n_env_vars)),
                      min=c(rep(0.05,n_extra),rep(-20,n_env_vars),rep(-10,n_env_vars)),
                      mean=c(rep(1,n_extra),rep(-15,n_env_vars),rep(-3,n_env_vars)),
                      sd=rep(1000,n_rows))

  return(pars_var)
}
#-------------------------------------------------------------------------------
#' @title pars_fixed_setup
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param sero_template TBA
#' @param case_template TBA
#' @param vacc_data TBA
#' @param pop_data TBA
#' @param years_data TBA
#' @param year0 TBA
#' @param time_inc TBA
#' @param mode_start TBA
#' @param start_SEIRV TBA
#' @param ... TBA
#'
#' @export
#'
pars_fixed_setup <- function(sero_template=list(),case_template=list(), vacc_data = list(),
                             pop_data = list(), years_data = c(), year0 = 1940,  time_inc = 1.0,
                             mode_start = 0, start_SEIRV = NULL, ...){

  vars_extra=list(...)
  assertthat::assert_that(all(names(vars_extra) %in% extra_param_names))

  assertthat::assert_that(mode_start %in% c(0, 1, 2), msg = "mode_start must have value 0, 1 or 2")
  if(mode_start==2){assertthat::assert_that(is.null(start_SEIRV)==FALSE)}
  assertthat::assert_that(length(pop_data[1, , 1]) > 1, msg = "Need population data for multiple years")
  assertthat::assert_that(length(pop_data[1, 1, ]) > 1, msg = "Need population data for multiple age groups")
  n_regions = length(pop_data[, 1, 1])
  n_years = length(pop_data[1, , 1]) - 1
  N_age = length(pop_data[1, 1, ])
  assertthat::assert_that(length(vacc_data[1, , 1]) == n_years + 1,
                          msg = "Population and vaccination data must be for same time periods")
  assertthat::assert_that(length(vacc_data[, 1, 1])==n_regions, msg = "No. regions in population and vaccination data must match")
  assertthat::assert_that(length(vacc_data[1, 1, ]) == N_age, msg="No. age groups in population and vaccination data must match")
  assertthat::assert_that(years_data[1] >= year0, msg = "First data year must be greater than or equal to year0")
  #Need 2 years of leeway - TBA
  assertthat::assert_that(max(years_data) + 2 - year0 <= n_years, msg = "Period of years_data must lie within population data with 2 years at end")
  assertthat::assert_that(time_inc %in% c(1, 2.5, 5), msg = "time_inc must have value 1, 2.5 or 5 days")
  inv_365 = 1.0/365.0
  pts_year = 365.0/time_inc
  n_t_pts = n_years*pts_year
  region=""

  #Initial SEIRV values
  P0 = array(0, dim = c(n_regions, N_age))
  for(i in 1:n_regions){for(j in 1:N_age){ P0[i,j] = max(1.0, pop_data[i,1, j]) }}

  if(mode_start == 2){ #Use supplied SEIRV data
    S_0 = start_SEIRV$S
    E_0 = start_SEIRV$E
    I_0 = start_SEIRV$I
    R_0 = start_SEIRV$R
    V_0 = start_SEIRV$V
  } else {
    E_0 = I_0 = array(0, dim = c(n_regions, N_age))
    vacc_initial = array(vacc_data[, 1, ],dim=c(n_regions,N_age))
    V_0 = P0*vacc_initial
    if(mode_start == 0){ #No initial immunity
      S_0 = P0*(1.0 - vacc_initial)
      R_0 = E_0
    } #if mode_start==1, S_0 and R_0 set based on FOI/R0 - see packer_setup
  }

  #Population/vaccination
  dP1_all = dP2_all = vacc_rates = array(NA, dim = c(n_regions, N_age, n_years))
  for(n_year in 1:n_years){
    for(i in 1:n_regions){
      for(j in 1:N_age){
        dP1_all[i,j,n_year] = max(1.0, pop_data[i, n_year + 1, j])*inv_365
        dP2_all[i,j,n_year] = max(1.0, pop_data[i, n_year, j])*inv_365
        if(j == 1){
          vacc_rates[i,j,n_year] = vacc_data[i, n_year + 1, j]*inv_365
        } else {
          vacc_rates[i,j,n_year] = max(0.0, vacc_data[i, n_year + 1, j] - vacc_data[i, n_year, j - 1])*inv_365
        }
      }
    }
  }

  #Region and output grouping
  sero_region_groups=sort(unique(sero_template$region))
  sero_regions=YEP::regions_breakdown(sero_region_groups)
  case_region_groups=sort(unique(case_template$region))
  case_regions=YEP::regions_breakdown(case_region_groups)
  regions_all=YEP::regions_breakdown(c(sero_template$region,case_template$region))
  n_regions=length(regions_all)
  i_sero_regions=i_case_regions=rep(0,n_regions)
  i_sero_regions[regions_all %in% sero_regions]=1
  i_case_regions[regions_all %in% case_regions]=1

  n_lines_sero=nrow(sero_template)
  unique_lines_sero=unique_lines_case=c(1)
  for(i in 2:n_lines_sero){
    unique=FALSE
    if(sero_template$region[i] %in% sero_template$region[unique_lines_sero]){
      if(sero_template$age_min[i] %in% subset(sero_template[c(1:(i-1)),],
                                              region==sero_template$region[i])$age_min == FALSE ||
         sero_template$age_max[i] %in% subset(sero_template[c(1:(i-1)),],
                                              region==sero_template$region[i])$age_max == FALSE){unique=TRUE}
    } else {
      unique=TRUE
    }
    if(unique){unique_lines_sero=append(unique_lines_sero,i,after=length(unique_lines_sero))}
  }
  n_sero_pts=length(unique_lines_sero)
  n_lines_case=nrow(case_template)
  for(i in 2:n_lines_case){
    unique=FALSE
    if(case_template$region[i] %in% case_template$region[unique_lines_case]==FALSE){unique=TRUE}
    if(unique){unique_lines_case=append(unique_lines_case,i,after=length(unique_lines_case))}
  }
  n_case_pts=length(unique_lines_case)
  region_index_sero=array(0,dim=c(n_sero_pts,n_regions))
  region_index_case=array(0,dim=c(n_case_pts,n_regions))
  sero_i_age_min=sero_i_age_max=sero_vc_factor=rep(NA,n_sero_pts)
  for(i in 1:n_sero_pts){
    regions_sero_pt=YEP::regions_breakdown(sero_template$region[unique_lines_sero[i]])
    for(j in 1:n_regions){
      if(regions_all[j] %in% regions_sero_pt){region_index_sero[i,j]=1}
    }
    sero_i_age_min[i]=sero_template$age_min[unique_lines_sero[i]]+1
    sero_i_age_max[i]=sero_template$age_max[unique_lines_sero[i]]+1
    sero_vc_factor[i]=sero_template$vc_factor[unique_lines_sero[i]]
  }
  for(i in 1:n_case_pts){
    regions_case_pt=YEP::regions_breakdown(case_template$region[unique_lines_case[i]])
    for(j in 1:n_regions){
      if(regions_all[j] %in% regions_case_pt){region_index_case[i,j]=1}
    }
  }
  output=list(n_r = n_regions, N_age = N_age, mode_start = mode_start,
              n_years = n_years, n_t_pts = n_t_pts, year0 = year0, time_inc = time_inc,
              t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious,
              E_0 = E_0, I_0 = I_0, V_0 = V_0, dP1_all = dP1_all, dP2_all = dP2_all,
              vacc_rate_daily = vacc_rates,
              n_sero_pts=n_sero_pts,n_case_pts=n_case_pts, sero_vc_factor=sero_vc_factor,
              sia_min=sero_i_age_min, sia_max=sero_i_age_max,
              region_index_sero = region_index_sero, region_index_case = region_index_case,
              sero_regions=i_sero_regions,case_regions=i_case_regions)
  if(mode_start %in% c(0,2)){
    output$S_0 = S_0
    output$R_0 = R_0
  }else{
    output$P0 = P0
  }
  for(name in extra_param_names){
    if(is.null(vars_extra[[name]])==FALSE){output[[name]]=vars_extra[[name]]}
  }

  return(output)
}
#-------------------------------------------------------------------------------
#' @title fit_data_setup
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param sero_template TBA
#' @param case_template TBA
#' @param year0 TBA
#' @param time_inc TBA
#' @param region_index_sero TBA
#' @param region_index_case TBA
#'
#' @export
#'
fit_data_setup <- function(sero_template = list(),case_template = list(), year0=1940, time_inc=5.0,
                           region_index_sero = list(), region_index_case = list()){

  year=region=age_min=age_max=0
  years_data <- sort(unique(c(sero_template$year,case_template$year)))
  i_year_begin=years_data[1] - year0 + 1
  i_year_end=max(years_data) + 1 - year0
  time_pts=(years_data-year0+1)*(365/time_inc)
  sero_region_groups=sort(unique(sero_template$region))
  case_region_groups=sort(unique(case_template$region))
  regions_all=YEP::regions_breakdown(c(sero_template$region,case_template$region))
  n_regions=length(regions_all)
  assertthat::assert_that(n_regions==dim(region_index_sero)[2])
  assertthat::assert_that(n_regions==dim(region_index_case)[2])
  n_sero_pts=dim(region_index_sero)[1]
  n_case_pts=dim(region_index_case)[1]

  n_t_pts=length(time_pts)
  sero_data_list1=sero_data_list2=case_data_list1=list()
  for(i in 1:n_t_pts){
    sero_data_list1[[i]]=sero_data_list2[[i]]=as.numeric(rep(NA,n_sero_pts))
    case_data_list1[[i]]=as.numeric(rep(NA,n_case_pts))
    year_pt=years_data[i]
    sero_subset=subset(sero_template,year==year_pt)
    if(nrow(sero_subset)>0){
      for(j in 1:n_sero_pts){
        regions=regions_all[which(region_index_sero[j,]==1)]
        region_group=paste(regions,collapse=",")
        sero_subset2=subset(sero_subset,region==region_group)
        sero_subset2=subset(sero_subset2,age_min==sero_template$age_min[j])
        sero_subset2=subset(sero_subset2,age_max==sero_template$age_max[j])
        sero_data_list1[[i]][j]=sero_subset2$positives[1]
        sero_data_list2[[i]][j]=sero_subset2$samples[1]
      }
    }
    case_subset=subset(case_template,year==year_pt)
    if(nrow(case_subset)>0){
      for(j in 1:n_case_pts){
        regions=regions_all[which(region_index_case[j,]==1)]
        region_group=paste(regions,collapse=",")
        case_data_list1[[i]][j]=case_subset$cases[case_subset$region==region_group]
      }
    }
  }
  fit_data=data.frame(time=time_pts,obs_sero_positives=I(sero_data_list1), obs_sero_samples=I(sero_data_list2),
                      obs_case_values=I(case_data_list1))

  return(fit_data)
}
#-------------------------------------------------------------------------------
# TODO - Make vector of scalar parameter names variable
#' @title packer_setup
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param pars_fixed TBA
#' @param env_covar_values TBA
#' @param mode_time TBA
#'
#' @export
#'
packer_setup <- function(pars_fixed = list(), env_covar_values = list(), mode_time=0){

  #TODO - add assert_that functions?
  assertthat::assert_that(length(dim(env_covar_values))==3)
  n_env_vars=dim(env_covar_values)[1]
  assertthat::assert_that(dim(env_covar_values)[1]==n_env_vars)
  assertthat::assert_that(mode_time %in% c(0:5), msg = "mode_time must be an integer between 0 and 5")
  assertthat::assert_that(all(env_covar_values >= 0.0))
  n_regions=pars_fixed$n_r
  assertthat::assert_that(dim(env_covar_values)[2]==n_regions)
  time_inc=pars_fixed$time_inc
  pts_year = 365.0/time_inc
  n_years=pars_fixed$n_years
  n_t_pts = n_years*pts_year
  inv_365 = 1.0/365.0
  n_req = switch(mode_time + 1, 1, n_years, 12, pts_year, n_years*12, n_t_pts)
  assertthat::assert_that(dim(env_covar_values)[3]==n_req)
  date_values = switch(mode_time + 1,
                       rep(1, n_t_pts),
                       sort(rep(c(1:n_years), pts_year)),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12),
                       1 + (floor(time_inc*c(0:(n_t_pts - 1))) %% pts_year),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12) +
                         (12*sort(rep(c(1:n_years) - 1,pts_year))),
                       c(1:n_t_pts))

  packer <- monty_packer(scalar = c("p_rep_severe","p_rep_death"),
                         array = list("log_FOI_coeffs" = n_env_vars,
                                      "log_R0_coeffs" = n_env_vars),
                         fixed = pars_fixed,
                         process = function(p){
                           FOI_spillover=colSums(exp(p$log_FOI_coeffs)*env_covar_values)
                           #if(m_FOI_BRA %in% p){} #TODO - add in m_FOI_BRA
                           R0=colSums(exp(p$log_R0_coeffs)*env_covar_values)
                           FOI_spillover_t=FOI_spillover[,date_values]
                           R0_t=R0[,date_values]
                           params_list=list(FOI_spillover=FOI_spillover_t,R0=R0_t)
                           if(pars_fixed$mode_start==1){
                             vacc_initial=pars_fixed$V_0
                             N_age=pars_fixed$N_age
                             P0=pars_fixed$P0
                             S_0 = R_0 = array(0, dim = c(n_regions, N_age))
                             ages = c(1:N_age) - 1
                             for(i in 1:n_regions){
                               R0_year0=mean(R0_t[i,c(1:pts_year)])
                               FOI_spillover_year0=mean(FOI_spillover_t[i,c(1:pts_year)])
                               if(R0_year0 <= 1.0){
                                 FOI_estimate = FOI_spillover_year0*365.0
                               } else {
                                 estimation_results = nlm(YEP::imm_fraction_function, p =  - 4, R0_year0, ages, P0[i,]/sum(P0[i,]))
                                 FOI_estimate = min(0.1, (FOI_spillover_year0*365.0) + exp(estimation_results$estimate))
                               }
                               herd_immunity = 1.0 - (exp( - FOI_estimate*(ages + 0.5)))

                               for(j in 1:N_age){
                                 if(vacc_initial[i,j]<herd_immunity[j]){
                                   R_0[i,j] = P0[i,j]*(herd_immunity[j] - vacc_initial[i,j])
                                   S_0[i,j] = P0[i,j]*(1.0 - herd_immunity[j])
                                 } else {
                                   R_0[i,j] = 0
                                   S_0[i,j] = P0[i,j]*(1.0 - vacc_initial[i,j])
                                 }
                               }
                             }
                             params_list$S_0=S_0
                             params_list$R_0=R_0
                           }
                           return(params_list)
                         })
  return(packer)
}
#-------------------------------------------------------------------------------
#' @title prior_setup
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param packer TBA
#' @param env_covar_values TBA
#' @param pars_var TBA
#' @param FOI_R0_prior_data TBA
#'
#' @export
#'
prior_setup <- function(packer=NULL, env_covar_values = list(), pars_var = list(), FOI_R0_prior_data = list()){

  #TODO - assert_that functions
  n_env_vars=dim(env_covar_values)[1]
  n_regions=dim(env_covar_values)[2]
  n_req=dim(env_covar_values)[3]
  cols1=which(grepl("log_FOI_coeffs",pars_var$name))
  cols2=which(grepl("log_R0_coeffs",pars_var$name))
  cols_add=which(c(1:nrow(pars_var)) %in% c(cols1,cols2) == FALSE)
  vars_extra=pars_var$name[cols_add]
  n_add=length(cols_add)
  env_covars_mean=array(NA,dim=c(n_env_vars,n_regions))
  for(i in 1:n_regions){
    env_covars_mean[,i]=rowMeans(array(env_covar_values[,i,],dim=c(n_env_vars,n_req)))
  }

  prior_function <- function(log_FOI_coeffs=rep(-15,5),log_R0_coeffs=rep(-5,5),
                             p_rep_severe=0,p_rep_death=0,p_severe_inf=0.12,p_death_severe_inf=0.39,
                             vaccine_efficacy=1.0,m_FOI_BRA=1.0,...){

    prior_lfc=prior_lrc=rep(0,n_env_vars)
    for(i in 1:n_env_vars){
      j=cols1[i]
      k=cols2[i]
      prior_lfc = log(truncdist::dtrunc(x=log_FOI_coeffs[i],spec="norm",a=pars_var$min[j],b=pars_var$max[j],
                                        mean=pars_var$mean[j],sd=pars_var$sd[j]))
      prior_lrc = log(truncdist::dtrunc(x=log_R0_coeffs[i],spec="norm",a=pars_var$min[k],b=pars_var$max[k],
                                        mean=pars_var$mean[k],sd=pars_var$sd[k]))
    }

    FOI_mean=R0_mean=rep(NA,n_regions)
    for(i in 1:n_regions){
      FOI_mean[i]=sum(env_covars_mean[,i]*exp(log_FOI_coeffs))
      R0_mean[i]=sum(env_covars_mean[,i]*exp(log_R0_coeffs))
    }
    prior_FOI=log(truncdist::dtrunc(x=FOI_mean,spec="norm",a=FOI_R0_prior_data$min[1],b=FOI_R0_prior_data$max[1],
                                    mean=FOI_R0_prior_data$mean[1],sd=FOI_R0_prior_data$sd[1]))
    prior_R0=log(truncdist::dtrunc(x=R0_mean,spec="norm",a=FOI_R0_prior_data$min[2],b=FOI_R0_prior_data$max[2],
                                   mean=FOI_R0_prior_data$mean[2],sd=FOI_R0_prior_data$sd[2]))

    prior_add=rep(0,n_add)
    for(i in 1:n_add){
      j=cols_add[i]
      prior_add[i]=log(truncdist::dtrunc(x=get(vars_extra[i]),spec="norm",a=pars_var$min[j],
                                         b=pars_var$max[j],mean=pars_var$mean[j],sd=pars_var$sd[j]))
    }

    return(sum(c(prior_lfc,prior_lrc,prior_FOI,prior_R0,prior_add)))
  }
  prior <- monty_model_function(prior_function,packer=packer)

  return(prior)
}
#-------------------------------------------------------------------------------
#' @title m_prelim_fit
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param fit_data TBA
#' @param packer TBA
#' @param prior TBA
#' @param FOI_R0_prior_data TBA
#' @param pars_var TBA
#' @param env_covar_values TBA
#' @param n_values TBA
#' @param n_particles TBA
#' @param n_threads TBA
#'
#' @export
#'
m_prelim_fit <- function(fit_data=list(), packer=NULL, prior=NULL, FOI_R0_prior_data = list(),
                         pars_var=list(), env_covar_values = list(),
                         n_values=10, n_particles=1,n_threads=1){

  filter <- dust_filter_create(generator = SEIRV_Model_mr04_fit, data = fit_data, time_start = 0,
                               n_particles = n_particles, n_threads = n_threads)
  likelihood <- dust_likelihood_monty(obj = filter, packer = packer)
  if(is.null(prior)){posterior <- likelihood}else{posterior<-likelihood+prior}

  n_env_vars=dim(env_covar_values)[1]
  n_regions=dim(env_covar_values)[2]
  n_req=dim(env_covar_values)[3]
  cols1=which(grepl("log_FOI_coeffs",pars_var$name))
  cols2=which(grepl("log_R0_coeffs",pars_var$name))

  param_sets=tgp::lhs(n_values,rect=matrix(c(pars_var$min,pars_var$max),ncol=2))
  explore=data.frame(array(NA,dim=c(n_values,dim(param_sets)[2]+1)))
  explore[,c(1:dim(param_sets)[2])]=param_sets
  colnames(explore)=c(paste0("param",c(1:dim(param_sets)[2])),"density")
  cat("\n")
  for(i in 1:n_values){
    time1=Sys.time()
    cat("\n",i)
    FOI_mean=R0_mean=rep(NA,n_regions)
    for(j in 1:n_regions){
      env_covars_mean=rowMeans(array(env_covar_values[,j,],dim=c(n_env_vars,n_req)))
      FOI_mean[j]=sum(env_covars_mean*exp(param_sets[i,cols1]))
      R0_mean[j]=sum(env_covars_mean*exp(param_sets[i,cols2]))
    }
    test1=any(FOI_mean>FOI_R0_prior_data$max[1])
    test2=any(R0_mean>FOI_R0_prior_data$max[2])
    if(any(test1,test2)){
      explore$density[i] <- -Inf
      cat("\n\t!!!")
    } else {
      explore$density[i] <- monty_model_density(posterior, parameters=param_sets[i,])
    }
    cat("\n",signif(param_sets[i,],4),"\n\t",explore$density[i])
    time2=Sys.time()
    cat("\n\t",time2-time1)
  }
  cat("\n")
  explore=explore[order(-explore$density),]

  return(explore)
}
#-------------------------------------------------------------------------------
#' @title m_sample
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param fit_data TBA
#' @param packer TBA
#' @param prior TBA
#' @param FOI_R0_prior_data TBA
#' @param pars_var TBA
#' @param initial TBA
#' @param v TBA
#' @param n_chains TBA
#' @param n_iterations TBA
#' @param n_particles TBA
#' @param n_threads TBA
#' @param parallel TBA
#' @param output_file TBA
#'
#' @export
#'
m_sample <- function(fit_data = list(), packer=NULL, prior=NULL, FOI_R0_prior_data = list(),
                     pars_var = list(), initial = list(), v = c(), n_chains = 1, n_iterations = 10,
                     n_particles = 1, n_threads = 1,  parallel = FALSE, output_file=""){

  #TODO - add assert_that functions
  n_params=nrow(pars_var)

  filter <- dust_filter_create(generator = SEIRV_Model_mr04_fit, data = fit_data, time_start = 0,
                               n_particles = n_particles, n_threads = n_threads)

  likelihood <- dust_likelihood_monty(obj = filter, packer = packer)
  if(is.null(prior)){posterior <- likelihood}else{posterior <- likelihood + prior}
  posterior$domain[,1]=pars_var$min
  posterior$domain[,2]=pars_var$max
  assertthat::assert_that(all(initial>=posterior$domain[,1]))
  assertthat::assert_that(all(initial<=posterior$domain[,2]))
  vcv <- matrix(c(rep(0,n_params^2)),n_params,n_params)
  for(i in 1:n_params){ vcv[i,i]=v[i] }
  sampler <- monty_sampler_random_walk(vcv = vcv, boundaries="reject",
                                       rerun_every=50, rerun_random=FALSE)
  #TODO - make rerun_every variable

  if(parallel){
    prior0=prior$density(initial[,1])#Seems to be necessary for some reason
    runner=monty_runner_callr(n_workers=n_chains)
  }else{runner=monty_runner_serial()}
  samples <- monty_sample(model = posterior,sampler = sampler,n_steps = n_iterations,
                          initial = initial, n_chains = n_chains, runner = runner)

  saveRDS(list(samples=samples,
               params=list(pars_var = pars_var, n_chains = n_chains,
                           n_iterations = n_iterations,
                           n_particles = n_particles, n_threads = n_threads)),
          file = output_file)

  return(samples)
}

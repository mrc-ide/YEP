# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YEP package - alternate versions
#-------------------------------------------------------------------------------
#' @title MCMC2
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param log_params_ini TBA
#' @param input_data TBA
#' @param obs_sero_data TBA
#' @param obs_case_data TBA
#' @param filename_prefix TBA
#' @param Niter TBA
#' @param type TBA
#' @param mode_start TBA
#' @param prior_settings TBA
#' @param dt TBA
#' @param n_reps TBA
#' @param enviro_data TBA
#' @param R0_fixed_values TBA
#' @param p_severe_inf TBA
#' @param p_death_severe_inf TBA
#' @param add_values TBA
#' @param deterministic TBA
#' @param mode_parallel TBA
#' @param cluster TBA
#' '
#' @export
#'
MCMC2 <- function(log_params_ini=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,filename_prefix="Chain",
                 Niter=1,type=NULL,mode_start=0,prior_settings=list(type="zero"),dt=1.0,n_reps=1,enviro_data=NULL,
                 R0_fixed_values=NULL,p_severe_inf=0.12,p_death_severe_inf=0.39,
                 add_values=list(vaccine_efficacy=1.0,p_rep_severe_af=1.0,p_rep_death_af=1.0,
                                 p_rep_severe_sa=1.0,p_rep_death_sa=1.0,m_FOI_Brazil=1.0),
                 deterministic=FALSE,mode_parallel="none",cluster=NULL){

  assert_that(is.logical(deterministic))
  assert_that(mode_start %in% c(0,1,3),msg="mode_start must have value 0, 1 or 3")
  n_params=length(log_params_ini)

  assert_that(all(names(add_values)==c("vaccine_efficacy","p_rep_severe_af","p_rep_death_af",
                                       "p_rep_severe_sa","p_rep_death_sa","m_FOI_Brazil"))) #Varied from default
  extra_estimated_params=c()
  for(var_name in names(add_values)){
    if(is.null(add_values[[var_name]])==TRUE){extra_estimated_params=append(extra_estimated_params,var_name)}
  }

  #Process input data to check that all regions with sero and/or case data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data (if used) and check that environmental data available for all regions
  input_data=input_data_process(input_data,obs_sero_data,obs_case_data)
  regions=names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions=length(regions)
  if(is.null(enviro_data)==FALSE){
    for(region in regions){assert_that(region %in% enviro_data$region)}
    enviro_data=subset(enviro_data,enviro_data$region %in% regions)
  }

  #Label parameters according to order and fitting type
  param_names=create_param_labels(type,input_data,enviro_data,extra_estimated_params)
  names(log_params_ini)=param_names

  #Run checks on inputs
  checks<-mcmc_checks(log_params_ini,n_regions,type,prior_settings,enviro_data,R0_fixed_values,add_values,
                      extra_estimated_params)
  if(prior_settings$type=="flat"){names(prior_settings$log_params_min)=names(prior_settings$log_params_max)=param_names}

  #Set up list of invariant parameter values to supply to other functions
  consts=list(type=type,mode_start=mode_start,prior_settings=prior_settings,dt=dt,n_reps=n_reps,enviro_data=enviro_data,
              R0_fixed_values=R0_fixed_values,p_severe_inf=p_severe_inf,p_death_severe_inf=p_death_severe_inf,
              add_values=add_values,extra_estimated_params=extra_estimated_params,
              deterministic=deterministic,mode_parallel=mode_parallel,cluster=cluster)

  #MCMC setup
  chain=chain_prop=posterior_current=posterior_prop=flag_accept=chain_cov_all=NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  log_params=log_params_ini
  chain_cov=1
  adapt=0
  posterior_value_current=-Inf

  #Iterative estimation
  for (iter in 1:Niter){

    #Propose new parameter values
    log_params_prop=param_prop_setup(log_params,chain_cov,adapt)

    #Calculate likelihood using single_posterior_calc function
    posterior_value_prop=single_posterior_calc2(log_params_prop,input_data,obs_sero_data,obs_case_data,consts)

    if(is.finite(posterior_value_prop)==FALSE) {
      p_accept = -Inf
    } else {
      p_accept = posterior_value_prop - posterior_value_current
      if(is.na(p_accept) ){ p_accept = -Inf}
    }

    ## accept/reject step:
    tmp = runif(1)
    if(tmp<min(exp(p_accept),1)) { # accept:
      log_params = log_params_prop
      posterior_value_current = posterior_value_prop
      accept = 1
    } else { # reject:
      accept = 0
    }

    #save current step
    chain = rbind(chain, log_params)
    chain_prop=rbind(chain_prop,log_params_prop)
    posterior_current=rbind(posterior_current,posterior_value_current)
    posterior_prop=rbind(posterior_prop,posterior_value_prop)
    flag_accept = rbind(flag_accept, accept)
    chain_cov_all = rbind(chain_cov_all,max(chain_cov))

    #Set output headings
    if(iter==1){
      colnames(chain)=colnames(chain_prop)=names(log_params_ini)
      for(i in 1:n_params){colnames(chain_prop)[i]=paste("Test_",colnames(chain_prop)[i],sep="")}
      colnames(posterior_current) = "posterior_current"
      colnames(posterior_prop) = "posterior_prop"
      colnames(flag_accept) = "flag_accept"
      colnames(chain_cov_all) = "chain_cov_all"
    }

    #Output chain to file every 10 iterations; start new file every 10,000 iterations
    if (iter %% 10 == 0){
      if (iter %% 10000 == 0){fileIndex  = iter/10000}

      if(fileIndex>=10){
        filename=paste0(filename_prefix,fileIndex,".csv")
      } else {
        filename=paste0(filename_prefix,"0",fileIndex,".csv")
      }
      if(file.exists(filename)==FALSE){file.create(filename)}
      lines=min((fileIndex * 10000+1),iter):iter

      data_out<-cbind(posterior_current,posterior_prop,exp(chain),flag_accept,exp(chain_prop),chain_cov_all)[lines,]
      if(fileAccess(filename,2)==0){write.csv(data_out,filename,row.names=FALSE)}
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }
  }

  #Get final parameter values
  param_out=exp(log_params)
  names(param_out)=names(log_params_ini)

  return(param_out)
}
#-------------------------------------------------------------------------------
#' @title single_posterior_calc2
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param log_params_prop TBA
#' @param input_data TBA
#' @param obs_sero_data TBA
#' @param obs_case_data TBA
#' @param consts TBA
#'
#' @export
#'
single_posterior_calc2 <- function(log_params_prop=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                                  consts=list()){

  regions=input_data$region_labels
  n_regions=length(regions)

  #Get additional values and calculate associated priors
  vaccine_efficacy=p_rep_severe_af=p_rep_death_af=p_rep_severe_sa=p_rep_death_sa=m_FOI_Brazil=1.0
  prior_add=0
  for(var_name in names(consts$add_values)){
    if(is.numeric(consts$add_values[[var_name]])==FALSE){
      i=match(var_name,names(log_params_prop))
      value=exp(as.numeric(log_params_prop[i]))
      assign(var_name,value)
      if(consts$prior_settings$type=="norm"){
        prior_add=prior_add+log(dtrunc(value,"norm",a=0,b=1,mean=consts$prior_settings$norm_params_mean[i],
                                       sd=consts$prior_settings$norm_params_sd[i]))
      } else {
        if(consts$prior_settings$type=="flat"){
          if(value<consts$prior_settings$log_params_min[i] || value>consts$prior_settings$log_params_max[i]){prior_add=-Inf}
        }
      }
    } else {
      assign(var_name,consts$add_values[[var_name]])
    }
  }

  #If additional values give finite prior, get FOI and R0 values and calculate associated prior
  if(is.finite(prior_add)){
    FOI_R0_data=mcmc_FOI_R0_setup(consts$type,consts$prior_settings,regions,log_params_prop,
                                  consts$enviro_data,consts$R0_fixed_values)
    FOI_values=FOI_R0_data$FOI_values
    for(n_region in 1:n_regions){
      if(substr(regions[n_region],1,3)=="BRA"){FOI_values[n_region]=FOI_values[n_region]*m_FOI_Brazil}
    }
    R0_values=FOI_R0_data$R0_values
    prior_prop=FOI_R0_data$prior+prior_add
  }else{prior_prop=-Inf}

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    if(is.null(obs_sero_data)){sero_like_values=0}
    if(is.null(obs_case_data)){cases_like_values=deaths_like_values=0}

    #Generate modelled data over all regions
    dataset <- Generate_Dataset2(input_data,FOI_values,R0_values,obs_sero_data,obs_case_data,vaccine_efficacy,consts$p_severe_inf,
                                consts$p_death_severe_inf,p_rep_severe_af,p_rep_death_af,p_rep_severe_sa,p_rep_death_sa,
                                consts$mode_start,start_SEIRV=NULL,consts$dt,
                                consts$n_reps,consts$deterministic,consts$mode_parallel,consts$cluster)

    #Likelihood of observing serological data
    if(is.null(obs_sero_data)==FALSE){
      sero_like_values=sero_data_compare(dataset$model_sero_values,obs_sero_data)
    }
    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data)==FALSE){
      cases_like_values=case_data_compare(dataset$model_case_values,obs_case_data$cases)
      if(is.null(obs_case_data$deaths)==FALSE){
        deaths_like_values=case_data_compare(dataset$model_death_values,obs_case_data$deaths)
      } else {deaths_like_values=0}
    }

    posterior=prior_prop+mean(c(sum(sero_like_values,na.rm=TRUE),sum(cases_like_values,na.rm=TRUE),
                                sum(deaths_like_values,na.rm=TRUE)),na.rm=TRUE)

  } else {posterior=-Inf}

  return(posterior)
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit2
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param n_iterations TBA
#' @param n_param_sets TBA
#' @param n_bounds TBA
#' @param type TBA
#' @param log_params_min TBA
#' @param log_params_max TBA
#' @param input_data TBA
#' @param obs_sero_data TBA
#' @param obs_case_data TBA
#' @param mode_start TBA
#' @param prior_settings TBA
#' @param dt TBA
#' @param n_reps TBA
#' @param enviro_data TBA
#' @param R0_fixed_values TBA
#' @param p_severe_inf TBA
#' @param p_death_severe_inf TBA
#' @param add_values TBA
#' @param deterministic TBA
#' @param mode_parallel TBA
#' @param cluster TBA
#' '
#' @export
#'
mcmc_prelim_fit2 <- function(n_iterations=1,n_param_sets=1,n_bounds=1,type=NULL,log_params_min=NULL,
                            log_params_max=NULL,input_data=list(),obs_sero_data=list(),obs_case_data=list(),
                            mode_start=0,prior_settings=list(type="zero"),dt=1.0,n_reps=1,enviro_data=NULL,R0_fixed_values=c(),
                            p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                            add_values=list(vaccine_efficacy=1.0,p_rep_severe_af=1.0,p_rep_death_af=1.0,
                                            p_rep_severe_sa=1.0,p_rep_death_sa=1.0,m_FOI_Brazil=1.0),
                            deterministic=TRUE,mode_parallel="none",cluster=NULL){

  #TODO - Add assertthat functions
  assert_that(mode_start %in% c(0,1,3),msg="mode_start must have value 0, 1 or 3")
  assert_that(length(log_params_min)==length(log_params_max),msg="Parameter limit vectors must have same lengths")
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(prior_settings$type %in% c("zero","exp","norm"))

  best_fit_results=list()
  n_params=length(log_params_min)

  #Get additional values
  extra_estimated_params=c()
  add_value_names=names(add_values)
  assert_that(add_value_names==c("vaccine_efficacy","p_rep_severe_af","p_rep_death_af",
                                 "p_rep_severe_sa","p_rep_death_sa","m_FOI_Brazil")) #Varied from default
  for(var_name in add_value_names){
    if(is.null(add_values[[var_name]])==TRUE){extra_estimated_params=append(extra_estimated_params,var_name)}
  }
  param_names=create_param_labels(type,input_data,enviro_data,extra_estimated_params)

  #TODO - Additional assert_that checks
  assert_that(length(param_names)==n_params)
  names(log_params_min)=names(log_params_max)=param_names
  xlabels=param_names
  for(i in 1:n_params){xlabels[i]=substr(xlabels[i],1,15)}
  ylabels=10^c(-8,-6,-4,-3,-2,-1,0,1)
  par(mar=c(6,2,1,1))
  ylim=c(min(log_params_min),max(log_params_max))

  for(iteration in 1:n_iterations){
    cat("\nIteration: ",iteration,"\n",sep="")
    all_param_sets <- lhs(n=n_param_sets,rect=cbind(log_params_min,log_params_max))
    results=data.frame()
    consts=list(type=type,mode_start=mode_start,prior_settings=prior_settings,
                dt=dt,n_reps=n_reps,enviro_data=enviro_data,R0_fixed_values=R0_fixed_values,
                p_severe_inf = p_severe_inf, p_death_severe_inf = p_death_severe_inf,add_values=add_values,
                deterministic=deterministic,mode_parallel=mode_parallel,cluster=cluster)

    for(set in 1:n_param_sets){
      cat("\n\tSet: ",set,sep="")
      log_params_prop=all_param_sets[set,]

      cat("\n\tParams: ",signif(log_params_prop,3))

      names(log_params_prop)=param_names
      posterior_value=single_posterior_calc(log_params_prop,input_data,obs_sero_data,obs_case_data,consts)
      results<-rbind(results,c(set,exp(log_params_prop),posterior_value))
      if(set==1){colnames(results)=c("set",param_names,"posterior")}

      cat("\n\tPosterior likelihood = ",posterior_value,sep="")

    }
    results<-results[order(results$posterior,decreasing=TRUE), ]
    best_fit_results[[iteration]]=results

    log_params_min_new=log_params_max_new=rep(0,n_params)
    for(i in 1:n_params){
      log_params_min_new[i]=min(log(results[c(1:n_bounds),i+1]))
      log_params_max_new[i]=max(log(results[c(1:n_bounds),i+1]))
    }
    names(log_params_min_new)=names(log_params_max_new)=param_names

    matplot(x=c(1:n_params),y=log(t(results[c(1:n_bounds),c(1:n_params)+1])),type="p",pch=16,col=1,
            xaxt="n",yaxt="n",xlab="",ylab="",ylim=ylim)
    axis(side=1,at=c(1:n_params),labels=xlabels,las=2,cex.axis=0.7)
    axis(side=2,at=log(ylabels),labels=ylabels)
    matplot(x=c(1:n_params),y=log_params_min,type="l",col=1,lty=2,add=TRUE)
    matplot(x=c(1:n_params),y=log_params_max,type="l",col=1,lty=2,add=TRUE)
    matplot(x=c(1:n_params),y=log_params_min_new,type="l",col=2,add=TRUE)
    matplot(x=c(1:n_params),y=log_params_max_new,type="l",col=2,add=TRUE)

    log_params_min=log_params_min_new
    log_params_max=log_params_max_new
  }

  return(best_fit_results)
}
#-------------------------------------------------------------------------------
#' @title Generate_Dataset2
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param input_data TBA
#' @param FOI_values TBA
#' @param R0_values TBA
#' @param sero_template TBA
#' @param case_template TBA
#' @param vaccine_efficacy TBA
#' @param p_severe_inf TBA
#' @param p_death_severe_inf TBA
#' @param p_rep_severe_af TBA
#' @param p_rep_death_af TBA
#' @param p_rep_severe_sa TBA
#' @param p_rep_death_sa TBA
#' @param mode_start TBA
#' @param start_SEIRV TBA
#' @param dt TBA
#' @param n_reps TBA
#' @param deterministic TBA
#' @param mode_parallel TBA
#' @param cluster TBA
#' @param output_frame TBA
#' '
#' @export
#'
Generate_Dataset2 <- function(input_data = list(),FOI_values = c(),R0_values = c(),sero_template = NULL,
                              case_template = NULL, vaccine_efficacy = 1.0,
                              p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                              p_rep_severe_af = 1.0,p_rep_death_af = 1.0,p_rep_severe_sa = 1.0,p_rep_death_sa = 1.0,
                              mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE,
                              mode_parallel = "none",cluster = NULL,output_frame=FALSE){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(any(is.null(sero_template)==FALSE,is.null(case_template)==FALSE),msg="Need at least one template")
  if(is.null(sero_template)==FALSE){
    assert_that(all(c("region","year","age_min","age_max","samples","vc_factor") %in% names(sero_template)))
  }
  if(is.null(case_template)==FALSE){
    assert_that(all(c("region","year") %in% names(case_template)))
    assert_that(p_severe_inf>=0.0 && p_severe_inf<=1.0,msg="Severe infection rate must be between 0-1")
    assert_that(p_death_severe_inf>=0.0 && p_death_severe_inf<=1.0,
                msg="Fatality rate of severe infections must be between 0-1")
    assert_that(all(c(p_rep_severe_af,p_rep_severe_sa,p_rep_death_af,p_rep_death_sa)>=0.0))
    assert_that(all(c(p_rep_severe_af,p_rep_severe_sa,p_rep_death_af,p_rep_death_sa)<=1.0))
  }
  assert_that(mode_parallel %in% c("none","pars_multi","clusterMap"))
  if(mode_parallel=="clusterMap"){assert_that(is.null(cluster)==FALSE)}

  #Prune input data based on regions
  regions=regions_breakdown(c(sero_template$region,case_template$region))
  input_data=input_data_truncate(input_data,regions)
  n_regions=length(input_data$region_labels)

  #Cross-reference templates with input regions
  if(is.null(sero_template)==FALSE){
    xref_sero=template_region_xref(sero_template,input_data$region_labels)
    sero_line_list=xref_sero$line_list
  } else {
    xref_sero=data.frame(year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))
    sero_line_list=rep(NA,n_regions)
  }
  if(is.null(case_template)==FALSE){
    country_list_af=c("AGO","BDI","BEN","BFA","CAF","CIV","CMR","COD","COG","ERI","ETH","GAB","GHA","GIN","GMB","GNB","GNQ",
                      "KEN","LBR","MLI","MRT","NER","NGA","RWA","SDN","SEN","SLE","SOM","SSD","TCD","TGO","TZA","UGA","ZMB")
    country_list_sa=c("ARG","BOL","BRA","COL","ECU","GUF","GUY","PER","PRY","SUR","VEN")
    countries=substr(regions,1,3)
    assert_that(all(countries %in% c(country_list_af,country_list_sa)))
    xref_case=template_region_xref(case_template,input_data$region_labels)
    case_line_list=xref_case$line_list
  } else {
    xref_case=data.frame(year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))
    case_line_list=rep(NA,n_regions)
  }
  year_data_begin=year_end=rep(NA,length(input_data$region_labels))
  for(i in 1:length(year_data_begin)){
    year_data_begin[i]=min(xref_sero$year_data_begin[i],xref_case$year_data_begin[i])
    year_end[i]=max(xref_sero$year_end[i],xref_case$year_end[i])
  }

  inv_reps=1/n_reps
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions to be modelled")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions to be modelled")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(sero_template)){model_sero_data=NULL} else {
    blank=rep(0,nrow(sero_template))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  }
  if(is.null(case_template)){model_case_values=model_death_values=NA} else {
    model_case_values=model_death_values=rep(0,nrow(case_template))
  }

  #Set up vector of output types to get from model if needed
  if(mode_parallel %in% c("none","clusterMap")){
    output_types=rep(NA,n_regions)
    for(n_region in 1:n_regions){
      if(is.na(case_line_list[[n_region]][1])==FALSE){
        if(is.na(sero_line_list[[n_region]][1])==FALSE){
          output_types[n_region]="case+sero"
        } else{
          output_types[n_region]="case"
        }
      } else {output_types[n_region]="sero"}
    }
  }

  #Model all regions in parallel if parallel modes in use
  if(mode_parallel=="pars_multi"){
    years_data_all=c(min(year_data_begin):max(year_end))
    if(is.null(sero_template)==FALSE){if(is.null(case_template)==FALSE){output_type="case+sero"} else {output_type="sero"}
    } else {output_type="case"}
    model_output_all=Model_Run_Multi_Input(FOI_spillover = FOI_values,R0 = R0_values,
                                           vacc_data = input_data$vacc_data, pop_data = input_data$pop_data,
                                           years_data = years_data_all, start_SEIRV=start_SEIRV,output_type = output_type,
                                           year0 = input_data$years_labels[1],mode_start = mode_start,
                                           vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                           n_threads = n_reps*n_regions,deterministic = deterministic)
  }
  if(mode_parallel=="clusterMap"){
    vacc_data_subsets=pop_data_subsets=years_data_sets=list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]]=input_data$vacc_data[n_region,,]
      pop_data_subsets[[n_region]]=input_data$pop_data[n_region,,]
      years_data_sets[[n_region]]=c(year_data_begin[n_region]:year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV=rep(NA,n_regions)}
    model_output_all=clusterMap(cl = cluster,fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                vacc_data = vacc_data_subsets,pop_data = pop_data_subsets,
                                years_data = years_data_sets, start_SEIRV = start_SEIRV, output_type = output_types,
                                MoreArgs=list(year0 = input_data$years_labels[1],mode_start = mode_start,
                                              vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,
                                              n_threads = 1 ,deterministic = deterministic))
  }
  #if(mode_parallel=="hybrid") #Potential future option combining parallelization types

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel=="none"){
      #cat("\n\t\tBeginning modelling region ",input_data$region_labels[n_region])
      model_output = Model_Run(FOI_spillover = FOI_values[n_region],R0 = R0_values[n_region],
                               vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                               years_data = c(year_data_begin[n_region]:year_end[n_region]),
                               start_SEIRV=start_SEIRV[[n_region]],output_type = output_types[n_region],
                               year0 = input_data$years_labels[1],mode_start = mode_start,
                               vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                               deterministic = deterministic)
      #cat("\n\t\tFinished modelling region ",n_region)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts=length(model_output$year)

    #Compile case data if needed
    if(is.na(case_line_list[[n_region]][1])==FALSE){
      case_line_list_region=case_line_list[[n_region]]
      years_case=case_template$year[case_line_list_region]
      n_lines=length(case_line_list_region)

      #Get reporting probabilities based on region
      if(countries[n_region] %in% country_list_af){
        p_rep_severe=p_rep_severe_af
        p_rep_death=p_rep_death_af
      } else {
        p_rep_severe=p_rep_severe_sa
        p_rep_death=p_rep_death_sa
      }

      for(n_rep in 1:n_reps){
        rep_cases=rep_deaths=rep(0,n_lines)
        for(n_line in 1:n_lines){
          #TODO - Set to adjust reporting probabilities based on region
          pts=c(1:t_pts)[model_output$year==years_case[n_line]]
          infs=sum(model_output$C[n_rep,pts])
          if(deterministic){
            severe_infs=floor(infs)*p_severe_inf
            deaths=severe_infs*p_death_severe_inf
            rep_deaths[n_line]=round(deaths*p_rep_death)
            rep_cases[n_line]=rep_deaths[n_line]+round((severe_infs-deaths)*p_rep_severe)

          } else {
            severe_infs=rbinom(1,floor(infs),p_severe_inf)
            deaths=rbinom(1,severe_infs,p_death_severe_inf)
            rep_deaths[n_line]=rbinom(1,deaths,p_rep_death)
            rep_cases[n_line]=rep_deaths[n_line]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
          }
        }

        model_case_values[case_line_list_region]=model_case_values[case_line_list_region]+rep_cases
        model_death_values[case_line_list_region]=model_death_values[case_line_list_region]+rep_deaths
      }
    }

    #Compile seroprevalence data if necessary
    if(is.na(sero_line_list[[n_region]][1])==FALSE){
      sero_line_list_region=sero_line_list[[n_region]]
      for(n_rep in 1:n_reps){
        sero_results=sero_calculate2(sero_template[sero_line_list_region,],model_output,n_rep)
        model_sero_data$samples[sero_line_list_region]=model_sero_data$samples[sero_line_list_region]+sero_results$samples
        model_sero_data$positives[sero_line_list_region]=model_sero_data$positives[sero_line_list_region]+sero_results$positives
      }
    }
  }

  if(is.null(sero_template)==FALSE){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}
  if(is.null(case_template)==FALSE){
    model_case_values=model_case_values*inv_reps
    model_death_values=model_death_values*inv_reps
  }

  if(output_frame) { #Output complete frames of data
    return(list(model_sero_data=data.frame(region=sero_template$region,year=sero_template$year,
                                           age_min=sero_template$age_min,age_max=sero_template$age_max,
                                           samples=sero_template$samples,positives=sero_template$samples*model_sero_data$sero),
                model_case_data=data.frame(region=case_template$region,year=case_template$year,
                                           cases=model_case_values,deaths=model_death_values)))
  } else { #Minimal output for MCMC
    return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
                model_death_values=model_death_values))
  }
}


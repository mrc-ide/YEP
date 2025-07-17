# Alternate version with cumulative annual calculation and run over multiple regions
# Version calculating both serological and case data and applying distribution
# FOI and R0 calculated in odin2 from environmental covariates and coefficients

#Parameters---------------------------------------------------------------------
time_inc <- parameter() #Time increment in days
n_r <- parameter() #number of regions
region_index_sero <- parameter() #Groupings of regions for which to output sero data
region_index_case <- parameter() #Groupings of regions for which to output case data
sero_regions <- parameter() #0/1 flag indicating which regions need serological data output
case_regions <- parameter() #0/1 flag indicating which regions need case data output
n_env_vars <- parameter() #number of environmental covariates
env_covar_values <- parameter() #Values of environmental covariates by variable, region and time point
log_FOI_coeffs <- parameter() #Log coefficients of environmental covariates used to calculate FOI_spillover by variable
log_R0_coeffs <- parameter() #Log coefficients of environmental covariates used to calculate R0 by variable
n_sero_pts <- parameter() #number of serology data points at each time point
n_case_pts <- parameter() #number of case data points at each time point
t_incubation <- parameter() #Length in days of yellow fever incubation period in mosquito vectors
t_latent <- parameter() #Length in days of latent period in humans exposed to yellow fever
t_infectious <- parameter() #Length of infectious period in humans with yellow fever
N_age <- parameter() #Number of age categories
vacc_rate_daily <- parameter() #Daily rate of vaccination by age and year
vaccine_efficacy <- parameter() #Proportion of vaccinations which successfully protect the recipient
sero_vc_factor <- parameter() #TBA
sia_min <- parameter() #TBA
sia_max <- parameter() #TBA
p_severe_inf <- parameter() #TBA
p_death_severe_inf <- parameter() #TBA
p_rep_severe <- parameter() #TBA
p_rep_death <- parameter() #TBA

#Initial conditions-------------------------------------------------------------
year0 <- parameter()  #Starting year
S_0 <- parameter() #Susceptible population by age group at start
E_0 <- parameter() #Exposed population by age group at start
I_0 <- parameter() #Infectious population by age group at start
R_0 <- parameter() #Recovered population by age group at start
V_0 <- parameter() #Vaccinated population by age group at start
dP1_all <- parameter() #Daily increase in number of people by age group (people arriving in group due to age etc.)
dP2_all <- parameter() #Daily decrease in number of people by age group (people leaving group due to age etc.)
n_years <- parameter() #Number of years for which model to be run
n_t_pts <- parameter() #Total number of time points
Pmin <- 1.0e-99 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
rate1 <- time_inc/(t_incubation+t_latent) # Rate of transfer from E to I
rate2 <- time_inc/t_infectious # Rate of transfer from I to R
p_rep <- p_severe_inf*((p_death_severe_inf*p_rep_death)+((1.0-p_death_severe_inf)*p_rep_severe)) #TBA

FOI_components[1:n_env_vars,1:n_r]=exp(log_FOI_coeffs[i])*env_covar_values[i,j,t_pt]
R0_components[1:n_env_vars,1:n_r]=exp(log_R0_coeffs[i])*env_covar_values[i,j,t_pt]
t_pt <- day/time_inc #Number of time points passed
FOI_spillover[1:n_r] <- sum(FOI_components[,i]) #Spillover force of infection (per day) at each time point
R0[1:n_r] <- sum(R0_components[,i]) #Basic reproduction number for human-human transmission at each time point
beta[1:n_r] <- (R0[i]*time_inc)/t_infectious #Daily exposure rate
FOI_sum[1:n_r] <-  min(FOI_max, beta[i]*(sum(I[i,])/P_tot[i]) + (FOI_spillover[i]*time_inc)) #Total force of infection
year_i <- floor(day/365)+1 #Number of years since start, as integer
flag_year <- if(as.integer(day+time_inc) %% 365 == 0) 1 else 0

dP1[1:n_r,1:N_age] <- dP1_all[i,j,year_i]*time_inc #Increase in population by age group over 1 time increment
dP2[1:n_r,1:N_age] <- dP2_all[i,j,year_i]*time_inc #Decrease in population by age group over 1 time increment
E_new[1:n_r,1:N_age] <- Binomial(as.integer(S[i,j]), FOI_sum[i]) #New exposed individuals by age group
I_new[1:n_r,1:N_age] <- E[i,j]*rate1     #New infectious individuals by age group
R_new[1:n_r,1:N_age] <- I[i,j]*rate2     #New recovered individuals by age group

P_nV[1:n_r,1:N_age] <- S[i,j] + R[i,j] #Total vaccine-targetable population by age group
inv_P_nV[1:n_r,1:N_age] <- 1.0/P_nV[i,j]
P[1:n_r,1:N_age] <- P_nV[i,j] + V[i,j] #Total population by age group (excluding E+I)
P_tot[1:n_r] <- sum(P[i, ]) #Total overall population (excluding E+I)
inv_P[1:n_r,1:N_age] <- 1.0/P[i,j]

vacc_rate[1:n_r,1:N_age] <- vacc_rate_daily[i,j,year_i]*vaccine_efficacy*time_inc*P[i,j] #Total no. vaccinations by age

#Updates to output values at each time increment--------------------------------
update(day) <- day + time_inc
update(year) <- year_i + year0 - 1
update(S[1:n_r,1]) <- max(Pmin, S[i,1] - E_new[i,1] - vacc_rate[i,1]*S[i,1]*inv_P_nV[i,1] + dP1[i,1] - (dP2[i,1]*S[i,1]*inv_P[i,1]))
update(S[1:n_r,2:N_age]) <- max(Pmin, S[i,j] - E_new[i,j] - vacc_rate[i,j]*S[i,j]*inv_P_nV[i,j] + (dP1[i,j]*S[i,j-1]*inv_P[i,j-1]) - (dP2[i,j]*S[i,j]*inv_P[i,j]))
update(E[1:n_r,1:N_age]) <- max(Pmin, E[i,j] + E_new[i,j] - I_new[i,j])
update(I[1:n_r,1:N_age]) <- max(Pmin, I[i,j] + I_new[i,j] - R_new[i,j])
update(R[1:n_r,1]) <- max(Pmin, R[i,1] + R_new[i,1] - vacc_rate[i,1]*R[i,1]*inv_P_nV[i,1] - (dP2[i,1]*R[i,1]*inv_P[i,1]))
update(R[1:n_r,2:N_age]) <- max(Pmin, R[i,j] + R_new[i,j] - vacc_rate[i,j]*R[i,j]*inv_P_nV[i,j] + (dP1[i,j]*R[i,j-1]*inv_P[i,j-1]) - (dP2[i,j]*R[i,j]*inv_P[i,j]))
update(V[1:n_r,1]) <- max(Pmin, V[i,1] + vacc_rate[i,1] - (dP2[i,1]*V[i,1]*inv_P[i,1]))
update(V[1:n_r,2:N_age]) <- max(Pmin, V[i,j] + vacc_rate[i,j] + (dP1[i,j]*V[i,j-1]*inv_P[i,j-1]) - (dP2[i,j]*V[i,j]*inv_P[i,j]))
update(R_cu[1:n_r,1:N_age]) <- if(sero_regions[i]==0) 0 else if(flag_year==1) 0 else 
  R_cu[i,j] + R[i,j]
update(SEIR_cu[1:n_r,1:N_age]) <- if(sero_regions[i]==0) 0 else if(flag_year==1) 0 else 
  SEIR_cu[i,j] + S[i,j]+E[i,j]+I[i,j]+R[i,j]
update(V_cu[1:n_r,1:N_age]) <- if(sero_regions[i]==0) 0 else if(flag_year==1) 0 else 
  V_cu[i,j] + V[i,j]
update(infs_cu[1:n_r,1:N_age]) <- if(case_regions[i]==0) 0 else if(flag_year==1) 0 else 
  infs_cu[i,j] + I_new[i,j]
update(R_an[1:n_sero_pts,1:n_r]) <- if(flag_year==0) 0 else if(region_index_sero[i,j]==0) 0 else 
  sum(R_cu[j,sia_min[i]:sia_max[i]]) + sum(R[j,sia_min[i]:sia_max[i]])
update(SEIR_an[1:n_sero_pts,1:n_r]) <- if(flag_year==0) 0 else if(region_index_sero[i,j]==0) 0 else 
  sum(SEIR_cu[j,sia_min[i]:sia_max[i]]) + sum(S[j,sia_min[i]:sia_max[i]])+sum(E[j,sia_min[i]:sia_max[i]])+sum(I[j,sia_min[i]:sia_max[i]])+sum(R[j,sia_min[i]:sia_max[i]])
update(V_an[1:n_sero_pts,1:n_r]) <- if(flag_year==0) 0 else if(region_index_sero[i,j]==0) 0 else 
  sum(V_cu[j,sia_min[i]:sia_max[i]]) + sum(V[j,sia_min[i]:sia_max[i]])
update(infs_an[1:n_case_pts,1:n_r]) <- if(flag_year==0) 0 else if(region_index_case[i,j]==0) 0 else 
  sum(infs_cu[j,1:N_age]) + sum(I_new[j,1:N_age])
update(output_sero[1:n_sero_pts]) <- if(sero_vc_factor[i]==0) sum(R_an[i,])/sum(SEIR_an[i,]) else 
  ((1.0-sero_vc_factor[i])*(sum(R_an[i,])/sum(SEIR_an[i,]))) +(sero_vc_factor[i]*((sum(R_an[i,])+sum(V_an[i,]))/(sum(SEIR_an[i,])+sum(V_an[i,]))))
update(output_case[1:n_case_pts]) <- Binomial(as.integer(sum(infs_an[i,])),p_rep) #TODO - calculate both cases and deaths

#Initial values of updated variables--------------------------------------------
initial(day) <- time_inc
initial(year) <- year0
initial(S[1:n_r,1:N_age]) <- S_0[i,j]
initial(E[1:n_r,1:N_age]) <- E_0[i,j]
initial(I[1:n_r,1:N_age]) <- I_0[i,j]
initial(R[1:n_r,1:N_age]) <- R_0[i,j]
initial(V[1:n_r,1:N_age]) <- V_0[i,j]
initial(R_cu[1:n_r,1:N_age]) <- 0
initial(SEIR_cu[1:n_r,1:N_age]) <- 0
initial(V_cu[1:n_r,1:N_age]) <- 0
initial(infs_cu[1:n_r,1:N_age]) <- 0
initial(R_an[1:n_sero_pts,1:n_r]) <- 0
initial(SEIR_an[1:n_sero_pts,1:n_r]) <- 0
initial(V_an[1:n_sero_pts,1:n_r]) <- 0
initial(infs_an[1:n_case_pts,1:n_r]) <- 0
initial(output_sero[1:n_sero_pts]) <- 0
initial(output_case[1:n_case_pts]) <- 0

#Dimensions---------------------------------------------------------------------
#Updated values
dim(S) <- c(n_r, N_age)
dim(E) <- c(n_r, N_age)
dim(I) <- c(n_r, N_age)
dim(R) <- c(n_r, N_age)
dim(V) <- c(n_r, N_age)
dim(R_cu) <- c(n_r, N_age)
dim(SEIR_cu) <- c(n_r, N_age)
dim(V_cu) <- c(n_r, N_age)
dim(infs_cu) <- c(n_r, N_age)
dim(R_an) <- c(n_sero_pts,n_r)
dim(SEIR_an) <- c(n_sero_pts,n_r)
dim(V_an) <- c(n_sero_pts,n_r)
dim(infs_an) <- c(n_case_pts,n_r)
dim(output_sero) <- n_sero_pts
dim(output_case) <- n_case_pts

#Calculated values
dim(FOI_components) <- c(n_env_vars,n_r)
dim(R0_components) <- c(n_env_vars,n_r)
dim(beta) <- n_r
dim(FOI_spillover) <- n_r
dim(R0) <- n_r
dim(FOI_sum) <- n_r
dim(dP1) <- c(n_r, N_age)
dim(dP2) <- c(n_r, N_age)
dim(E_new) <- c(n_r, N_age)
dim(I_new) <- c(n_r, N_age)
dim(R_new) <- c(n_r, N_age)
dim(P_nV) <- c(n_r, N_age)
dim(inv_P_nV) <- c(n_r, N_age)
dim(P) <- c(n_r, N_age)
dim(P_tot) <- n_r
dim(inv_P) <- c(n_r, N_age)
dim(vacc_rate) <- c(n_r, N_age)

#Inputs
dim(region_index_sero) <- c(n_sero_pts,n_r)
dim(region_index_case) <- c(n_case_pts,n_r)
dim(sero_regions) <- n_r
dim(case_regions) <- n_r
dim(env_covar_values) <- c(n_env_vars,n_r, n_t_pts)
dim(log_FOI_coeffs) <- n_env_vars
dim(log_R0_coeffs) <- n_env_vars
dim(S_0) <- c(n_r, N_age)
dim(E_0) <- c(n_r, N_age)
dim(I_0) <- c(n_r, N_age)
dim(R_0) <- c(n_r, N_age)
dim(V_0) <- c(n_r, N_age)
dim(dP1_all) <- c(n_r, N_age, n_years)
dim(dP2_all) <- c(n_r, N_age, n_years)
dim(vacc_rate_daily) <- c(n_r, N_age, n_years)
dim(sero_vc_factor) <- n_sero_pts
dim(sia_min) <- n_sero_pts
dim(sia_max) <- n_sero_pts

#Distribution-------------------------------------------------------------------
#obs_sero_values <- data()
#dim(obs_sero_values) <- n_sero_pts
#obs_sero_values[] ~ Poisson(output_sero[i])
obs_sero_positives <- data()
obs_sero_samples <- data()
dim(obs_sero_positives) <- n_sero_pts
dim(obs_sero_samples) <- n_sero_pts
obs_sero_positives[] ~ Binomial(size = obs_sero_samples[i], prob = output_sero[i])
obs_case_values <- data()
dim(obs_case_values) <- n_case_pts
#obs_case_values[] ~ Poisson(output_case[i])
case_population <- data()
dim(case_population) <- n_case_pts
obs_case_values[] ~ NegativeBinomial(size = case_population[i], mu = output_case[i])

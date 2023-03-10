# SEIRV (Susceptible, Exposed, Infectious, Recovered, Vaccinated) yellow fever model,
# incorporating the force of infection of spillover from sylvatic/non-human primate reservoirs (which can also
# represent case importation) and the reproduction number for human-human transmission. Returns SEIRV data at each
# time point (separated by increment dt) and also numbers of new infections and total force of infection at each
# time point.

#Parameters
dt <- user() #Time increment in days
FOI_spillover <- user() #Spillover force of infection (per day)
R0 <- user() #Basic reproduction number
N_age <- user() #Number of age categories
vacc_rate_annual[,] <- user() #Daily rate of vaccination by age and year
vaccine_efficacy <- user() #Proportion of vaccinations which successfully protect the recipient

#initial conditions
year0 <- user()  #Starting year
Sus0[] <- user() #Susceptible population by age group at start
Exp0[] <- user() #Exposed population by age group at start
Inf0[] <- user() #Infectious population by age group at start
Rec0[] <- user() #Recovered population by age group at start
Vac0[] <- user() #Vaccinated population by age group at start
Cas0[] <- user() #Daily cases population by age group at start
dP1_all[,] <- user() #Daily increase in number of people by age group (people arriving in group due to age etc.)
dP2_all[,] <- user() #Daily decrease in number of people by age group (people leaving group due to age etc.)
n_years <- user() #Number of years for which model to be run

Pmin <- 0 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
I_lag[1:N_age] <- delay(I[i], t_incubation)
beta <- (R0*dt)/t_infectious #Daily exposure rate
FOI_sum <-  min(FOI_max, beta*(sum(I_lag)/P_tot) + (FOI_spillover*dt)) #Total force of infection
year_i <- floor(((step+1)*dt)/365) + 1 #Number of years since start, as integer
dP1[1:N_age] <- dP1_all[i, as.integer(year_i)]*dt #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, as.integer(year_i)]*dt #Decrease in population by age group over 1 time increment

E_new[1:N_age] <- S[i]*FOI_sum #New exposed individuals by age group (deterministic)
I_new[1:N_age] <- delay(E_new[i], t_latent)     #New infectious individuals by age group
R_new[1:N_age] <- delay(I_new[i], t_infectious)     #New recovered individuals by age group
P_nV[1:N_age] <- S[i] + R[i] #Total vaccine-targetable population by age group
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- P_nV[i] + V[i] #Total population by age group
P_tot <- sum(P) #Total overall population
inv_P[1:N_age] <- 1.0/P[i]
vacc_rate[1:N_age] <- vacc_rate_annual[i,as.integer(year_i)]*vaccine_efficacy*dt*P[i] #Total no. vaccinations by age

#Updates to output values at each time increment
update(day) <- day + dt
update(year) <- year_i + year0 - 1
update(FOI_total) <- FOI_sum
update(S[1]) <- max(Pmin,S[1] - E_new[1] - vacc_rate[1]*S[1]*inv_P_nV[1] + dP1[1] - (dP2[1]*S[1]*inv_P[1]))
update(S[2:N_age]) <- max(Pmin,S[i] - E_new[i] - vacc_rate[i]*S[i]*inv_P_nV[i] + (dP1[i]*S[i-1]*inv_P[i-1]) - (dP2[i]*S[i]*inv_P[i]))
update(E[1:N_age]) <- max(Pmin,E[i] + E_new[i] - I_new[i])
update(I[1:N_age]) <- max(Pmin,I[i] + I_new[i] - R_new[i])
update(R[1]) <- max(Pmin,R[1] + R_new[1] - vacc_rate[1]*R[1]*inv_P_nV[1] - (dP2[1]*R[1]*inv_P[1]))
update(R[2:N_age]) <- max(Pmin,R[i] + R_new[i] - vacc_rate[i]*R[i]*inv_P_nV[i] + (dP1[i]*R[i-1]*inv_P[i-1]) - (dP2[i]*R[i]*inv_P[i]))
update(V[1]) <- max(Pmin,V[1] + vacc_rate[1] - (dP2[1]*V[1]*inv_P[1]))
update(V[2:N_age]) <- max(Pmin,V[i] + vacc_rate[i] + (dP1[i]*V[i-1]*inv_P[i-1]) - (dP2[i]*V[i]*inv_P[i]))
update(C[1:N_age]) <- I_new[i]

#Initial values
initial(day) <- 0
initial(year) <- year0-1
initial(FOI_total) <- FOI_spillover
initial(S[1:N_age]) <- Sus0[i]
initial(E[1:N_age]) <- Exp0[i]
initial(I[1:N_age]) <- Inf0[i]
initial(R[1:N_age]) <- Rec0[i]
initial(V[1:N_age]) <- Vac0[i]
initial(C[1:N_age]) <- Cas0[i]

#Dimensions
dim(S) <- N_age
dim(E) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(V) <- N_age
dim(C) <- N_age

dim(I_lag) <- N_age
dim(dP1)<-N_age
dim(dP2)<-N_age
dim(E_new) <- N_age
dim(I_new) <- N_age
dim(R_new) <- N_age
dim(P_nV) <- N_age
dim(inv_P_nV) <- N_age
dim(P) <- N_age
dim(inv_P) <- N_age
dim(vacc_rate) <- N_age

dim(Sus0) <- N_age
dim(Exp0) <- N_age
dim(Inf0) <- N_age
dim(Rec0) <- N_age
dim(Vac0) <- N_age
dim(Cas0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_annual) <- c(N_age, n_years)

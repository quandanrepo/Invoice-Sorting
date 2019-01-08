#Read in data
data <- read.csv("data/1002_GB_IvatExtract.csv")

#load in library
library(tidyverse)

#Remove columns that are NULL and NA
# Can check indivual columns using all(is.na()), say all(is.na(data$VATCode))
# We want all columns that have NA values

col_remove_na <- which(sapply(data, function(x)all(is.na(x))))

#Subset the data, by removing the above columns
#use data_reduce just to keep original
data_reduce <- subset(data,select = -c(col_remove_na))

#Check columns with all NULL

col_remove_null <- which(sapply(data_reduce, function(x)all(x=="NULL")))

#Subset data
data_reduce <- subset(data_reduce,select = -c(col_remove_null))

#Check for columns with only 0 values

col_remove_zero <- which(sapply(data_reduce, function(x)all(x==0)))

#Subset data
data_reduce <- subset(data_reduce,select = -c(col_remove_zero))

#Still have many columns with NULL entires but have only a few entries, therefore not removed

#Example 1 - looking at invoice amounts over time
invoice_amounts <- data_reduce %>% 
  dplyr::mutate(date = as.Date(Transaction.Date)) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(invoice_amount = sum(Total.Value))

#Converting to data.frame to use with ggplot
df_invoice_amo <- as.data.frame(invoice_amounts)
#Simple ggplot
invoice_plot1 <- ggplot(data = df_invoice_amo,aes(x=date,y=invoice_amount))+geom_line(color="red")+geom_point()

#Example 2 - looking at the number of invoices over time
invoice_volumes <- data_reduce %>% 
  dplyr::mutate(date = as.Date(Transaction.Date)) %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(invoice_volume = n())

#Converting to data.frame to use with ggplot
df_invoice_vol <- as.data.frame(invoice_volumes)
#Simple ggplot
invoice_plot2 <- ggplot(data = df_invoice_vol,aes(x=date,y=invoice_volume))+geom_line(color="red")+geom_point()

#Arrange diagrams together
gridExtra::grid.arrange(invoice_plot1,invoice_plot2)

#Basic observations:
#Invoice amount peaks at similar times throughout the years which can be explained
#by the number of invoices (invoice volumes) handed in at the corresponding times
#Can be looked as a seaonality issue, have a closer look

#Have a look at the total invoices per month
invoice_vol_month <- data_reduce %>% 
  dplyr::mutate(transaction_month=month(as.Date(Transaction.Date))) %>% 
  dplyr::group_by(transaction_month) %>% 
  dplyr::summarise(invoice_volume =n())

#Check if this can be classified by a normal of student-t distribution

library(LaplacesDemon) #for student-t distribution

# Set seed
set.seed(1)

# Input data is invoice_vol_month
# Transform as monthly log returns
y <- diff(log(as.data.frame(invoice_vol_month)[,2]))

# Set priors
# Assume equal chance of being in either model i.e. Prior set to be equal -> p(k=1)=p(k=2)=1/2
model_prior <- c(1/2,1/2)
mu_prior <- c(0,1) #prior parameters for mu
sigma_prior <- c(1) #prior parameters for sigma #Assume that sigma is distribution by half-normal distribution
df_prior <- c(3,30) #prior parameters for df

# Need to write a function to evaluate half-normal distribution
dhalfnorm <- function(y,scale_halfnorm,log=c(TRUE,FALSE)){
  if(log=="TRUE"){
    den_halfnorm <- log(2*dnorm(y,0,scale_halfnorm))
  } else {
    den_halfnorm <- 2*dnorm(y,0,scale_halfnorm)
  }
  return(den_halfnorm)
}

# proposal parameters
proposal_df <- c(3,30) #hyperparameter choices for df

# Begin sampler in k=1 i.e. normal model
sims <- 100000 #Iterations
mu_list <- rep(NA,sims) #Keep mu information
sigma_list <- rep(NA,sims) #Keep sigma information
df_list <- rep(NA,sims) #Keep df information

#Proposalsd for Metropolis Hastings
mu_proposal_mh <- c(0,1)
sigma_proposal_mh <- c(0,1)
df_proposal_mh <- c(0,1)

model_list <- rep(NA,sims) #Keep model information
para <- c(0,1) #Initialise mu and sigma
model <- 1 #Start in Normal model
for (s in 1:sims) {
  if(s%%10000==0){
    print(s)
  }
  if(model==1){#We propose jump to Student-t
    new_para <- c(para[1],para[2],runif(1,proposal_df[1],proposal_df[2]))
    
    old_like <- sum(dnorm(y,para[1],para[2],log = TRUE))
    new_like <- sum(dst(y,new_para[1],new_para[2],new_para[3],log=TRUE))
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dhalfnorm(para[2],sigma_prior,log=TRUE)  
    new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(new_para[3],df_prior[1],df_prior[2],log = TRUE) +dhalfnorm(new_para[2],sigma_prior,log = TRUE)
    
    old_to_new <- dunif(new_para[3],proposal_df[1],proposal_df[2],log = TRUE)
    new_to_old <- log(1) #To reverse jump, we have a deterministic proposal
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)){
      para <- new_para
      model <- 2
    }
  }
  else{#We propose jump to normal
    new_para <- c(para[1],para[2])
    
    old_like <- sum(dst(y,para[1],para[2],para[3],log=TRUE))
    new_like <- sum(dnorm(y,new_para[1],new_para[2],log = TRUE))
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(para[3],df_prior[1],df_prior[2],log = TRUE) +dhalfnorm(para[2],sigma_prior,log = TRUE)
    new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) +dhalfnorm(new_para[2],sigma_prior,log = TRUE)
    
    old_to_new <- log(1) #To reverse jump, we have a deterministic proposal
    new_to_old <- dunif(para[3],proposal_df[1],proposal_df[2],log = TRUE)
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior+new_to_old-old_to_new)){
      para <- new_para
      model <- 1
      
    }
  }
  #resample parameters for either model
  #metropolis hasting sampling
  if(model==1){
    new_para <- c(para[1]+rnorm(1,mu_proposal_mh[1],mu_proposal_mh[2]),para[2])
    
    old_like <- sum(dnorm(y,para[1],para[2],log = TRUE))
    new_like <- sum(dnorm(y,new_para[1],new_para[2],log = TRUE))
    
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dhalfnorm(para[2],sigma_prior,log = TRUE)
    new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) + dhalfnorm(new_para[2],sigma_prior,log = TRUE)
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior)){
      para <- new_para
    }
  }
  else {
    new_para <- c(para[1]+rnorm(1,mu_proposal_mh[1],mu_proposal_mh[2]),para[2],para[3]+rnorm(1,df_proposal_mh[1],df_proposal_mh[2]))
    
    if(new_para[3]<3){
      new_prior <- -Inf
      new_like <- -Inf
    }
    else{
      new_prior <- dnorm(new_para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(new_para[3],df_prior[1],df_prior[2],log = TRUE) + dhalfnorm(new_para[2],sigma_prior,log = TRUE)
      new_like <- sum(dst(y,new_para[1],new_para[2],new_para[3],log = TRUE))
    }
    
    old_like <- sum(dst(y,para[1],para[2],para[3],log = TRUE))
    
    old_prior <- dnorm(para[1],mu_prior[1],mu_prior[2],log = TRUE) + dunif(para[3],df_prior[1],df_prior[2],log = TRUE) + dhalfnorm(para[2],sigma_prior,log = TRUE)
    
    if(runif(1)<exp(new_like+new_prior-old_like-old_prior)){
      para <- new_para
    }
  }
  
  model_list[s] <- model
  mu_list[s] <- para[1]
  sigma_list[s] <- para[2]
  df_list[s] <- para[3]
}

model_list_fac <- factor(model_list)
levels(model_list_fac) <- c("Normal","Student-t")
par(mfrow=c(2,1))
plot(model_list_fac)
hist(df_list,breaks = 300)

############Likely that it is Normal Data than Student-t
#Apply my changepoint model to it just to see if there is a change in invoice numbers:
#If there is no change, that implies that there is regularity in this period.

par(mfrow=c(2,2)) #have 4 plots in same space
plot(y,col=1,type = 'o')

# prior for parameter
mu_prior <- c(0,0.3) #Normal distribution #Mu is conditional on sigma i.e sigma^2/kappa #unless we choose unknown mean and known sigma (1)
sigma_prior <- c(5,0.6) #Use Scaled-Inverse-Chi-Squared for Conjugacy (df_para,scale_para) - mean of c(5,0.6) = 1
tau_prior <- function(tau,y,log=c(TRUE,FALSE)){ #Just assume uniform for now
  if (tau<1|tau>(length(y)-1)){
    return(0)
  } else {
    if (log) {
      return(log(1/(length(y)-1)))
    } else {
      return(1/(length(y)-1))
    }
  }
}

# proposal_sd_rj <- 0.3 #bigger jump for the reversible jump
# proposal_sd_wm <- 0.1 #within model sample of parameters
tau_lambda <- c(1,length(y)/2) # lambda for poisson jump 
large_jump <- 0.01 #prob of doing a large jump
sims <- 50000 
approx_samples <- 100 #Number of samples for monte carlo estimate
perc_burnin <- 0.1
thin_every <- 5 # Add thinning as there might be autocorrelation since alot of proposals will get rejected due to out of bounds proposals for tau
proposal_tau <- function(tau,tau_lambda,large_jump){
  if (runif(1)<large_jump) {
    return(tau+(((-1)^(sample(0:1,1)))*rpois(1,tau_lambda[2])))
  } else{
    return(tau+(((-1)^(sample(0:1,1)))*rpois(1,tau_lambda[1])))
  }
  # return(sample(1:(length(y)-1),1))
}

# Set up parameters to store
num_of_cps_list <- rep(NA,sims)
tau_list <- list()
segment_1 <- list() #to hold both mu and sigma
segment_2 <- list() #to hold both mu and sigma

# Initialisation
tau <- sample(1:(length(y)-1),1)
mu <- mean(y)
sigma <- sqrt(var(y))

#Unknown Mu, Unknown Sigma^2
# Do analytical solution first
marg_like=rep(NA,(length(y)-1))

# Before anything, define function for product of normals
prod_norm <- function(mu,sigma,y){prod(dnorm(y,mu,sigma))}

# Define function to generate and evaluate Scaled-Inverse-Chi-Squared variables, parameters (scale,degree_of_freedom)
rscaleinvchi <- function(n,df_para,scale_para){
  return((df_para*scale_para)/rchisq(n,df_para))
}
dscaleinvchi <- function(x,df_para,scale_para,log=c(TRUE,FALSE)){
  if (log) {
    return(log((scale_para*df_para/2)^(df_para/2)/(gamma(df_para/2))*exp((-df_para*scale_para)/(2*x))/(x^(1+(df_para/2)))))
  } else {
    return((scale_para*df_para/2)^(df_para/2)/(gamma(df_para/2))*exp((-df_para*scale_para)/(2*x))/(x^(1+(df_para/2))))
  } 
}

# Write function for Marginal likelihood
marg_like_mu_sigma <- function(y,mu_prior,sigma_prior){
  df_n <- sigma_prior[1]+length(y)
  kappa_n <- mu_prior[2]+length(y)
  mu_n <- (mu_prior[2]*mu_prior[1]+length(y)*mean(y))/kappa_n
  scale_n <- (1/df_n)*(sigma_prior[1]*sigma_prior[2]+sum((y-mean(y))^2)+((length(y)*mu_prior[2])/(mu_prior[2]+length(y)))*(mu_prior[1]-mean(y))^2)
  return((1/(pi^(length(y)/2)))*(sqrt(mu_prior[2]/kappa_n))*(gamma(df_n/2)/gamma(sigma_prior[1]/2))*((sigma_prior[1]*sigma_prior[2])^(sigma_prior[1]/2))/(df_n*scale_n)^(df_n/2))
}

# Do analytical solution first
marg_like=rep(NA,(length(y)-1))

for (i in 1:(length(y)-1)) {
  marg_like[i] <- marg_like_mu_sigma(y[1:i],mu_prior,sigma_prior)*marg_like_mu_sigma(y[(i+1):(length(y))],mu_prior,sigma_prior)
}
# marg_like <- marg_like[-1]
marg_like <- marg_like/sum(marg_like)

plot(marg_like)
lines(marg_like)

results <- data.frame("Analytical"=marg_like)


# Approximation Monte Carlo -----------------------------------------------
# Will be alot tricker to do a approximation since we need to do
# 1. Generate approx_sample sigma, and for each sigma, generate approx sample mu
# 2. Calculate the integration of both segments by taking mean

# Marginal Likelihood P(Y|M_1)
test_sigma_M1_1 <- sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2]))
test_mu_M1_1 <- list()
for (i in 1:length(test_sigma_M1_1)) {
  test_mu_M1_1[[i]] <- rnorm(approx_samples,mu_prior[1],sqrt(test_sigma_M1_1[i]^2/mu_prior[2]))
}
test_sigma_M1_2 <- sqrt(rscaleinvchi(approx_samples,sigma_prior[1],sigma_prior[2]))
test_mu_M1_2 <- list()
for (i in 1:length(test_sigma_M1_2)) {
  test_mu_M1_2[[i]] <- rnorm(approx_samples,mu_prior[1],sqrt(test_sigma_M1_2[i]^2/mu_prior[2]))
}
pM1 <- list()
for (j in 1:(length(y)-1)) { #number of changepoint models
  temp_approx <- rep(NA,approx_samples) #initialise the storage
  for (i in 1:approx_samples) { # integrate over both parameters
    temp_approx[i] <- log(mean(sapply(test_mu_M1_1[[i]],prod_norm,test_sigma_M1_1[i],y[1:j]))*mean(sapply(test_mu_M1_2[[i]],prod_norm,test_sigma_M1_2[i],y[(j+1):(length(y))])))
  }
  pM1[j] <- mean(exp(temp_approx)) #need to take exp
}
approx_posterior <- unlist(pM1)
# Normalise
# total <- pM0+pM1
# pM0 <- pM0/total #Not accurate due to the number of samples required for a good approximation
approx_posterior <- approx_posterior/sum(approx_posterior)#

plot(approx_posterior)
lines(approx_posterior)
# Store results in data frame
results <- cbind(results,"Approx"=approx_posterior)


# Gibbs Sampler + MH on tau ------------------------------------------------------

sigma_posterior_dist <- function(y,mu_prior,sigma_prior){
  df_n <- sigma_prior[1]+length(y)
  scale2_n <- (1/df_n)*(sigma_prior[1]*sigma_prior[2]+sum((y-mean(y))^2)+((length(y)*mu_prior[2])/(mu_prior[2]+length(y)))*(mu_prior[1]-mean(y))^2)
  return(rscaleinvchi(1,df_n,scale2_n))
}

mu_posterior_dist <- function(y,mu_prior,sigma){
  kappa_n <- mu_prior[2]+length(y)
  mu_n <- (mu_prior[2]*mu_prior[1]+length(y)*mean(y))/kappa_n
  return(rnorm(1,mu_n,sqrt((sigma^2)/(kappa_n))))
}

for (s in 1:sims) {
  if (s%%10000==0) {
    print(s)
  }
  #Resample parameters in segment 1 using gibbs sampler
  #sample sigma2 first, then mu conditional on sigma2
  temp_sigma2 <- sigma_posterior_dist(y[1:tau],mu_prior,sigma_prior)
  temp_mu <- mu_posterior_dist(y[1:tau],mu_prior,sqrt(temp_sigma2))
  segment_1[[s]] <- c(temp_mu,sqrt(temp_sigma2))
  
  #Resample parameters in segment 2 using gibbs sampler
  #sample sigma2 first, then mu conditional on sigma2
  temp_sigma2 <- sigma_posterior_dist(y[(tau+1):length(y)],mu_prior,sigma_prior)
  temp_mu <- mu_posterior_dist(y[(tau+1):length(y)],mu_prior,sqrt(temp_sigma2))
  segment_2[[s]] <- c(temp_mu,sqrt(temp_sigma2))
  
  #Resample tau using Metropolis hasting
  #Propose a new tau usings a symmetric proposal so no transition ratio
  pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
  
  while (pro_tau<1|pro_tau>(length(y)-1)) {
    pro_tau <- proposal_tau(tau,tau_lambda,large_jump)
  }
  
  old_like <- sum(dnorm(y[1:tau],segment_1[[s]][1],segment_1[[s]][2],log = TRUE)) + sum(dnorm(y[(tau+1):length(y)],segment_2[[s]][1],segment_2[[s]][2],log = TRUE))
  new_like <- sum(dnorm(y[1:pro_tau],segment_1[[s]][1],segment_1[[s]][2],log = TRUE)) + sum(dnorm(y[(pro_tau+1):length(y)],segment_2[[s]][1],segment_2[[s]][2],log = TRUE))
  
  old_prior <- tau_prior(tau,y,log = TRUE)
  new_prior <- tau_prior(pro_tau,y,log = TRUE)
  
  if (runif(1)<exp(new_like+new_prior-old_like-old_prior)) {
    tau <- pro_tau
  }
  
  tau_list[[s]] <- tau
}
tau_list <- tau_list[-(1:(perc_burnin*sims))]
tau_freq <- as.data.frame(table(factor(unlist(tau_list[seq(1,length(tau_list),thin_every)]))))
tau_den <- tau_freq[,2]/sum(tau_freq[,2])

# plot(marg_like)
# lines(marg_like)
plot(tau_den)
lines(tau_den)

results <- cbind(results,"Gibbs+MH"=c(tau_den))

### The results show that there are statistically no changepoint points in this monthly summary of invoice number
### i.e. there is no regime change for the number of invoices recieved  
## (the four plots are just analytically derived, monte carlo estimate, gibbs sampler(most flexible as can be adapted to non-parameteric dirichlet models for clusters))

##Not sure if this was any use:
##Look at the difference in transaction date and process date
##using select(data_reduce,contains("date")), we find the processing date column
process_date <- data_reduce %>% 
  dplyr:: select(Transaction.Date,ProcessingDate) %>% 
  dplyr:: mutate(t_date=as.Date(Transaction.Date),p_date=as.Date(ProcessingDate))
plot(process_date$p_date-process_date$t_date,type = 'l')

difference_dates <- process_date$p_date-process_date$t_date
minimum_process_to_transaction_time <- which(difference_dates<quantile(difference_dates,1/20))
hist(minimum_process_to_transaction_time)
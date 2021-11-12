# Joanna Harton
# Last Updateed: August 2021

library(survival)
library(NPP)
library(optmatch)


###### data generation functions #######
make_surv_sim <- function(trt.prop, trial.size, ehr.fraction, 
                          trt.effect, on.trial.effect=1, 
                          trial.cens.rate=0.1, ehr.cens.rate=0.4,
                          confounding.strength=c("none", "mild", "moderate", "strong")) {
  
  # Make the trial data
  trial <- make_trial_surv_sim(trt.prop, trial.size, 
                               trt.effect, on.trial.effect, 
                               trial.cens.rate,
                               confounding.strength)
  # Make the EHR data
  ehr <- make_ehr_surv_sim(trt.prop, trial.size, ehr.fraction,
                           trt.effect, on.trial.effect,
                           ehr.cens.rate,
                           confounding.strength)
  
  dat <- rbind(trial, ehr)
  dat$ID <- 1:nrow(dat)
  # Puts ID column first
  dat <- dat[,c(ncol(dat),1:ncol(dat)-1)]
  
  return(dat)
  
} 


make_trial_surv_sim <- function(trt.prop, trial.size, 
                                trt.effect, on.trial.effect, 
                                trial.cens.rate,
                                confounding.strength=c("none", "mild", "moderate", "strong")) {
  
  nsubj <- trial.size
  
  # Calculates number of treated and control subjects
  num.trt <- trt.prop*trial.size
  num.control <- trial.size-num.trt
  # Randomizes the order of the treated and control subjects in the dataset
  poss.vals <- c(rep(1, num.trt), rep(0, num.control))
  trt <- sample(poss.vals, size=trial.size, replace = FALSE)
  on.trial <- rep(1, trial.size)
  # Sets up empty vectors for covariates
  x1 <- rep(NA, nsubj)
  x2 <- rep(NA, nsubj)
  x3 <- rep(NA, nsubj)
  x4 <- rep(NA, nsubj)
  
  # Fills in values for covariates
  for (ii in 1:nsubj) {
    x1[ii] <- rbinom(1, 1, 0.5)
    x2[ii] <- rbinom(1, 1, 0.6)
    x3[ii] <- rnorm(1, mean=60, sd=5)-60
    x4[ii] <- rnorm(1, mean=21, sd=2)-21
  }
  
  # x1 is binomial (male) - males have higher hazard of death than females
  # x2 is binomial (college degree) - college educated have lower hazard of death than non-college educated (proxy for SES, etc.)
  # x3 is continuous (HDL cholesterol) - higher HDL has lower hazard of death than lower HDL
  # x4 is continuous (BMI) - higher BMI has higher hazard of death than lower BMI
  
  # Sets strength of association between covariates and outcome
  if (confounding.strength=="none" | confounding.strength=="None") {
    beta1 <- 1
    beta2 <- 1
    beta3 <- 1
    beta4 <- 1
  } else if (confounding.strength=="mild" | confounding.strength=="Mild") {
    beta1 <- 1.25
    beta2 <- 0.67
    beta3 <- 0.98
    beta4 <- 1.06
  } else if (confounding.strength=="moderate" | confounding.strength=="Moderate") {
    beta1 <- 1.75
    beta2 <- 0.5
    beta3 <- 0.96
    beta4 <- 1.12
  } else if (confounding.strength=="strong" | confounding.strength=="Strong") {
    beta1 <- 2.25
    beta2 <- 0.4
    beta3 <- 0.93
    beta4 <- 1.21
  } else {
    time("Choose one of the following for confounding strength: none, mild, moderate, strong")
  }
  
  # Generate the randomness for individuals
  u <- runif(nsubj, 0, 1)
  
  # Create the failure time from the treatment and on-trial effects as well as the 4 covariates
  failure.time <- -log(u)/exp(log(trt.effect)*trt+log(on.trial.effect)*on.trial+log(beta1)*x1+log(beta2)*x2+log(beta3)*x3+log(beta4)*x4)
  
  # Censoring time is dependent on being on-trial or in the real world
  censoring.time <- rep(NA, nsubj)
  for (ii in 1:length(censoring.time)) {
    censoring.time[ii] <- rexp(1, rate=trial.cens.rate)
  }
  
  event <- rep(0, nsubj)
  time <- rep(NA, nsubj)
  
  # Create an event indicator and a 'time' variable that is the min of cens and failure
  for (ii in 1:nsubj) {
    time[ii] <- min(failure.time[ii], censoring.time[ii])
    if (failure.time[ii]<=censoring.time[ii]) {
      event[ii] <- 1
    }
  }
  
  # Create dataframe with all variables
  dat <- data.frame(trt, on.trial, x1, x2, x3, x4, time, event)
  
  return(dat)
  
} 


make_ehr_surv_sim <- function(trt.prop, trial.size, ehr.fraction, 
                              trt.effect, on.trial.effect, 
                              ehr.cens.rate,
                              confounding.strength=c("none", "mild", "moderate", "strong")) {
  
  # Calculate how many EHR subjects we have available
  nsubj <- ehr.fraction*trial.size
  
  trt <- rep(0, nsubj)
  on.trial <- rep(0, nsubj)
  
  # Sets up empty vectors for covariates
  x1 <- rep(NA, nsubj)
  x2 <- rep(NA, nsubj)
  x3 <- rep(NA, nsubj)
  x4 <- rep(NA, nsubj)
  
  # Fills in values for covariates
  for (ii in 1:nsubj) {
    x1[ii] <- rbinom(1, 1, 0.55)
    x2[ii] <- rbinom(1, 1, 0.4)
    x3[ii] <- rnorm(1, mean=60, sd=10)-60
    x4[ii] <- rnorm(1, mean=23, sd=2)-21
  }
  
  
  # x1 is binomial (male) - males have higher hazard of death than females
  # x2 is binomial (college degree) - college educated have lower hazard of death than non-college educated (proxy for SES, etc.)
  # x3 is continuous (HDL cholesterol) - higher HDL has lower hazard of death than lower HDL
  # x4 is continuous (BMI) - higher BMI has higher hazard of death than lower BMI
  
  # Sets strength of association between covariates and outcome
  if (confounding.strength=="none" | confounding.strength=="None") {
    beta1 <- 1
    beta2 <- 1
    beta3 <- 1
    beta4 <- 1
  } else if (confounding.strength=="mild" | confounding.strength=="Mild") {
    beta1 <- 1.25
    beta2 <- 0.67
    beta3 <- 0.98
    beta4 <- 1.06
  } else if (confounding.strength=="moderate" | confounding.strength=="Moderate") {
    beta1 <- 1.75
    beta2 <- 0.5
    beta3 <- 0.96
    beta4 <- 1.12
  } else if (confounding.strength=="strong" | confounding.strength=="Strong") {
    beta1 <- 2.25
    beta2 <- 0.4
    beta3 <- 0.93
    beta4 <- 1.21
  } else {
    time("Choose one of the following for confounding strength: none, mild, moderate, strong")
  }
  
  # Generate the randomness for individuals
  u <- runif(nsubj, 0, 1)
  
  
  # Create the failure time from the treatment and on-trial effects as well as the 4 covariates
  failure.time <- -log(u)/exp(log(trt.effect)*trt+log(on.trial.effect)*on.trial+log(beta1)*x1+log(beta2)*x2+log(beta3)*x3+log(beta4)*x4)
  
  # Censoring time is dependent on being on-trial or in the real world
  censoring.time <- rep(NA, nsubj)
  for (ii in 1:length(censoring.time)) {
    censoring.time[ii] <- rexp(1, rate=ehr.cens.rate)
  }
  
  event <- rep(0, nsubj)
  time <- rep(NA, nsubj)
  
  # Create an event indicator and a 'time' variable that is the min of cens and failure
  for (ii in 1:nsubj) {
    time[ii] <- min(failure.time[ii], censoring.time[ii])
    if (failure.time[ii]<=censoring.time[ii]) {
      event[ii] <- 1
    }
  }
  
  # Create dataframe with all variables
  dat <- data.frame(trt, on.trial, x1, x2, x3, x4, time, event)
  
  
  return(dat)
  
} 





####### individual method functions #####

fit_ignore_external <- function(dat) {
  # only uses trial subjects
  dat_trial <- dat[dat$on.trial==1, ]

  cox_ie <- coxph(Surv(time, event)~trt, dat=dat_trial)
  
  cox.ie <- summary(cox_ie)
  
  my.estimates <- unname(c(est=cox.ie$coefficients[1]))
  my.ESS <- nrow(dat_trial)
  my.var <- cox.ie$coefficients[3]^2
  
  return(c(est=my.estimates, var=my.var, ESS=my.ESS, alpha=0))
}

fit_full_pooling <- function(dat) {
  # uses all trial and all EHR subjects
  cox_fp <- coxph(Surv(time, event)~trt, dat=dat)
  
  cox.fp <- summary(cox_fp)
  
  my.estimates <- unname(c(est=cox.fp$coefficients[1]))
  my.ESS <- nrow(dat)
  my.var <- cox.fp$coefficients[3]^2
  
  return(c(est=my.estimates, var=my.var, ESS=my.ESS, alpha=1))
}

fit_power_prior <- function(dat, alpha) {
  dat_trial <- dat[dat$on.trial==1, ]
  dat_ehr <- dat[dat$on.trial==0, ]
  
  # Give all trial subjects weight of 1
  dat_trial$weight <- 1
  # Give all EHR subjects weight of alpha
  dat_ehr$weight <- alpha
  
  dat <- rbind(dat_trial, dat_ehr)
  
  cox_pp <- coxph(Surv(time, event)~trt, dat=dat, weights = dat$weight)
  
  cox.pp <- summary(cox_pp)
  
  my.estimates <- unname(c(est=cox.pp$coefficients[1]))
  my.ESS <- sum(dat$weight)
  my.var <- cox.pp$coefficients[3]^2
  

  return(c(est=my.estimates, var=my.var, ESS=my.ESS, alpha=alpha))
}

fit_normalized_power_prior <- function(dat) {
  
  dat_trial <- dat[dat$on.trial==1 & dat$trt==0, ]
  dat_ehr <- dat[dat$on.trial==0, ]
  
  # Prep for estimating alpha hat
  time_trial <- dat_trial$time
  time_ext <- dat_ehr$time
  X_trial <- as.matrix(dat_trial[ , c("x1", "x2", "x3", "x4")])
  X_ext <- as.matrix(dat_ehr[ , c("x1", "x2", "x3", "x4")])
  
  # Estimate alpha hat
  out <- LMNPP_MCMC(y.Cur = time_trial, y.Hist = time_ext, x.Cur = X_trial, x.Hist = X_ext,
                    prior = list(a = 1.5, b = 0, mu0 = c(0, 0),
                                 Rinv = diag(100, nrow = 2),
                                 delta.alpha = 1, delta.beta = 1), MCMCmethod = 'RW',
                    rw.logit.delta = 2, nsample = 5000,
                    control.mcmc = list(delta.ini = NULL,
                                        burnin = 2000, thin = 5))
  
  
  dat_trial <- dat[dat$on.trial==1, ]
  dat_ehr <- dat[dat$on.trial==0, ]
  
  # Give all trial subjects weight of 1
  dat_trial$weight <- 1
  # Give all EHR subjects weight of alpha hat
  dat_ehr$weight <- mean(out$delta)
  
  dat <- rbind(dat_trial, dat_ehr)
  
  cox_npp <- coxph(Surv(time, event)~trt, dat=dat, weights = dat$weight)
  
  cox.npp <- summary(cox_npp)
  
  my.estimates <- unname(c(est=cox.npp$coefficients[1]))
  my.ESS <- sum(dat$weight)
  my.var <- cox.npp$coefficients[3]^2
  
  return(c(est=my.estimates, var=my.var, ESS=my.ESS, alpha=mean(out$delta)))
}

est_propensity_score <- function(dat) {
  
  # Estimate the propensity score with a main effects logistic model (data generating model)
  ps.model <- glm(on.trial~x1+x2+x3+x4, data=dat, family=binomial())
  dat$ps <- fitted(ps.model)
  
  for (ii in 1:nrow(dat)) {
    if (dat$ps[ii]>1) {
      dat$ps[ii] <- max(dat$ps[dat$ps<1])
    }
  }
  
  # Calculate IOW for all subjects
  dat$iow.weight <- 1
  for (ii in 1:nrow(dat)) {
    if (dat$on.trial[ii]==0) {
      dat$iow.weight[ii] <- dat$ps[ii]/(1-dat$ps[ii])
    }
  }

  trial.active <- dat[dat$on.trial==1 & dat$trt==1, ]
  trial.control <- dat[dat$on.trial==1 & dat$trt==0, ]
  external <- dat[dat$on.trial==0, ]
  
  dat2 <- rbind(trial.active, trial.control, external)

  return(dat2)
}

fit_LIN_method <- function(dat) {
  
  all.dat <- dat
  # Estimate propensity score via logistic regression with all trial patient (treated and control) and external controls
  trial.treated <- dat[dat$trt==1, ]
  trial.control <- dat[dat$on.trial==1 & dat$trt==0, ]
  external.control <- dat[dat$on.trial==0, ]
  
  # Match external controls with treated trial subjects so you have a 1:1 treated : control ratio
  dat.trt.external <- rbind(trial.treated, external.control)
  
  # Pair match using optimal matching (optmatch package) and then randomly select (n_T - n_C) external controls
  my.x <- dat.trt.external$ps
  names(my.x) <- rownames(dat.trt.external)
  
  my.z <- dat.trt.external$on.trial
  names(my.z) <- rownames(dat.trt.external)
  
  pair.matches <- pairmatch(x=my.x, z=my.z, controls=1, data=dat.trt.external)
  dat.trt.external$matches <- as.numeric(pair.matches)
  
  # Remove the unmatchables
  unmatchable <- dat.trt.external[is.na(dat.trt.external$matches)==TRUE, ]
  matchable <- dat.trt.external[is.na(dat.trt.external$matches)==FALSE, ]
  matchable.external <- matchable[matchable$on.trial==0, ]
  rownames(matchable.external) <- 1:nrow(matchable.external)
  
  # Determine how many external controls we need
  num.ext.needed <- nrow(dat[dat$trt==1, ])-nrow(dat[dat$on.trial==1 & dat$trt==0, ])
  
  chosen <- sample(1:nrow(matchable.external), num.ext.needed, replace = FALSE)
  selected.external <- matchable.external[chosen, ]
  selected.external$matches <- NULL
  selected.external.IDs <- as.data.frame(selected.external$ID)
  colnames(selected.external.IDs) <- "ID"
  selected.external.IDs$Lin.selected <- 1
  
  new.dat.pair.matched <- rbind(trial.treated, trial.control, selected.external)
  
  # Weight the external controls by their propensity score as the alpha parameter
  dat_trial <- new.dat.pair.matched[new.dat.pair.matched$on.trial==1, ]
  dat_ehr <- new.dat.pair.matched[new.dat.pair.matched$on.trial==0, ]
  
  dat_trial$weight <- 1
  dat_ehr$weight <- dat_ehr$ps
  
  dat <- rbind(dat_trial, dat_ehr)
  
  
  test <- merge(all.dat, selected.external.IDs, by="ID", all=TRUE)
  for (ii in 1:nrow(test)) {
    if (test$on.trial[ii]==0 & is.na(test$Lin.selected[ii])==TRUE) {
      test$Lin.selected[ii] <- 0
    }
  }
  
  cox_lin <- coxph(Surv(time, event)~trt, dat=dat, weights = dat$weight)
  
  cox.lin <- summary(cox_lin)
  
  my.estimates <- unname(c(est=cox.lin$coefficients[1]))
  my.ESS <- sum(dat$weight)
  my.var <- cox.lin$coefficients[3]^2
  
  return(list(out=c(est=my.estimates, var=my.var, ESS=my.ESS, alpha=NA), dat=test))
}

fit_data_adaptive_weighting <- function(dat) {
  all.dat <- dat
  
  dat_ehr <- dat[dat$on.trial==0, ]

  dat_trial <- dat[dat$on.trial==1, ]
  
  # Sort the EHR subjects in decreasing order
  new_dat_ehr <- dat_ehr[order(dat_ehr$ps, decreasing=T), ]
  # Determine how many EHR subjects are needed to make augmented trial 1:1
  num_ehr_needed <- nrow(dat_trial[dat_trial$trt==1,])-nrow(dat_trial[dat_trial$trt==0,])
  # Take number needed with largest propensity scores
  selected_ehr <- new_dat_ehr[1:num_ehr_needed, ]
  selected_ehr.IDs <- as.data.frame(selected_ehr$ID)
  colnames(selected_ehr.IDs) <- "ID"
  selected_ehr.IDs$DAW.selected <- 1
  
  # Scale the IOW so they sum to the number of EHR subjects selected
  selected_ehr$scaled.iow.weight <- (selected_ehr$iow.weight/sum(selected_ehr$iow.weight))*length(selected_ehr$iow.weight)
  dat_trial$scaled.iow.weight <- dat_trial$iow.weight
  
  dat <- rbind(dat_trial, selected_ehr)
  
  
  test <- merge(all.dat, selected_ehr.IDs, by="ID", all=TRUE)
  for (ii in 1:nrow(test)) {
    if (test$on.trial[ii]==0 & is.na(test$DAW.selected[ii])==TRUE) {
      test$DAW.selected[ii] <- 0
    }
  }
  ## 
  cox_daw <- coxph(Surv(time, event)~trt, dat=dat, weights = dat$scaled.iow.weight)
    
  cox.daw <- summary(cox_daw)
    
  my.estimates <- unname(c(est=cox.daw$coefficients[1]))
  my.ESS <- sum(dat$scaled.iow.weight)
  my.var <- cox.daw$coefficients[3]^2
  
  return(list(out=c(est=my.estimates, var=my.var, ESS=my.ESS, alpha=NA), dat=test))
}


###### Wrapper function  ####

# set up dataframe for results
empty.dat <- as.data.frame(matrix(NA, nrow=8, ncol=4))
colnames(empty.dat) <- c("est", "var", "ESS", "alpha")
rownames(empty.dat) <- c("ie", "fp", "pp_0.25", "pp_0.5", "pp_0.75", "npp", "lin", "daw")

# Set up nested list for dataframe to live in
out1 <- list(results=empty.dat, data=NA)

confound <- list(no.confounding=out1, mild.confounding=out1, strong.confounding=out1)
trt.effect <- list(trt.0.5=confound, trt.0.75=confound, trt.0.875=confound, trt.1=confound)
ehr.fraction <- list(ehr.frac.1=trt.effect, ehr.frac.5=trt.effect)
trial.size <- list(trial.100=ehr.fraction, trial.1000=ehr.fraction)
trt.prop <- list(trt.prop.0.67=trial.size, trt.prop.0.75=trial.size)

ALL <- trt.prop

all_methods <- function(ii) {
  
  # Possible values for parameters
  poss.trt.prop <- c(0.67, 0.75)
  poss.trial.size <- c(100, 1000)
  poss.ehr.fraction <- c(1, 5)
  poss.trt.effect <- c(0.5, 0.75, 0.875, 1)
  poss.confounding.strength <- c("none", "mild", "strong")
  
  # Nested loop running through all combinations
  for (aa in 1:length(poss.trt.prop)) {
    for (bb in 1:length(poss.trial.size)) {
      for (cc in 1:length(poss.ehr.fraction)) {
        for (dd in 1:length(poss.trt.effect)) {
          for (ff in 1:length(poss.confounding.strength)) {
            
            # Generate the data
            dat <- make_surv_sim(trt.prop=poss.trt.prop[aa], 
                                 trial.size=poss.trial.size[bb], 
                                 ehr.fraction=poss.ehr.fraction[cc], 
                                 trt.effect=poss.trt.effect[dd], 
                                 confounding.strength=poss.confounding.strength[ff])
            
            # Estimate the propensity score
            dat.w.ps <- est_propensity_score(dat)

            # Use each method
            ie <- tryCatch({fit_ignore_external(dat.w.ps)}, error = function(e) {rep(NA,4)})
            fp <- tryCatch({fit_full_pooling(dat.w.ps)}, error = function(e) {rep(NA,4)})
            pp_0.25 <- tryCatch({fit_power_prior(dat.w.ps, alpha=0.25)}, error = function(e) {rep(NA,4)})
            pp_0.5 <- tryCatch({fit_power_prior(dat.w.ps, alpha=0.5)}, error = function(e) {rep(NA,4)})
            pp_0.75 <- tryCatch({fit_power_prior(dat.w.ps, alpha=0.75)}, error = function(e) {rep(NA,4)})
            npp <- tryCatch({fit_normalized_power_prior(dat.w.ps)}, error = function(e) {rep(NA,4)})
            all.lin <- tryCatch({fit_LIN_method(dat.w.ps)}, error = function(e) {rep(NA,4)})
            lin <- all.lin$out
            all.daw <- tryCatch({fit_data_adaptive_weighting(all.lin$dat)}, error = function(e) {rep(NA,4)})
            daw <- all.daw$out
            
            out <- rbind(ie, fp, pp_0.25, pp_0.5, pp_0.75, npp, lin, daw)
            
            ALL[[aa]][[bb]][[cc]][[dd]][[ff]] <- list(results=out, data=all.daw$dat)
          }
        }
      }
    }
  }
  
  saveRDS(ALL, file=paste0("June_All_Methods_FREQ_COX/Replicate", ii, ".RDS"))
  
}


### Run ####

oldw <- getOption("warn")
options(warn=-1)

set.seed(8888)
# Randomize seed order
poss.seeds <- 1:100000
seeds <- sample(poss.seeds, size=100000, replace=FALSE)


overall <- function(i, wd = '/home/jograce/Simulations_DAW') {
  # set seed for each rep for reproducibility
  set.seed(seeds[i])
  all_methods(i)
  
}

results = parallel::mclapply(1:5000, overall, wd = '/home/jograce/Simulations_DAW', mc.cores = 50)

options(warn=oldw)


# Joanna Harton
# Last updated: August 2021

nsims <- 5000

# Estimated true marginal HR by running trial only method with a sample size of 1,000,000
marginal_HR <- readRDS("true_marginal_HRs.RDS")

# Function to make dataframe for all reps and methods together
gen_empty <- function() {
  
  empty.dat <- as.data.frame(matrix(NA, nrow=nsims, ncol=8))
  
  colnames(empty.dat) <- c("ie", "fp",
                           "pp_0.25", "pp_0.5", "pp_0.75",
                           "npp", "lin", "daw")
  
  confound <- list(no.confounding=empty.dat, mild.confounding=empty.dat, strong.confounding=empty.dat)
  trt.effect <- list(trt.0.5=confound, trt.0.75=confound, trt.0.875=confound, trt.1=confound)
  ehr.fraction <- list(ehr.frac.1=trt.effect, ehr.frac.5=trt.effect)
  trial.size <- list(trial.100=ehr.fraction, trial.1000=ehr.fraction)
  trt.prop <- list(trt.prop.0.67=trial.size, trt.prop.0.75=trial.size)
  ALL <- trt.prop
  
  
  return(ALL)
}

# Determines whether CI covers truth or not
obtain_ci_coverage <- function(lb, ub, my.num) {
  
  coverage <- rep(NA, length(lb))
  for (jj in 1:length(lb)) {
    if (is.na(lb[jj])==F & is.na(ub[jj])==F) {
      if ((lb[jj] <= my.num & ub[jj] >= my.num) | (ub[jj] <= my.num & lb[jj] >= my.num)) {coverage[jj] <- 1}
      else {coverage[jj] <- 0}
    }
  }
  
  return(coverage)
}


# Createe empty dataframes
my.log.hrs <- gen_empty()
my.hrs <- gen_empty()
my.vars <- gen_empty()
my.ESSs <- gen_empty()
my.coverage.one <- gen_empty()
my.coverage.truth <- gen_empty()
my.coverage.one.rounded <- gen_empty()
my.coverage.truth.rounded <- gen_empty()
my.alphas <- gen_empty()

# Possible values
poss.trt.prop <- c(0.67, 0.75)
poss.trial.size <- c(100, 1000)
poss.ehr.fraction <- c(1, 5)
poss.trt.effect <- c(0.5, 0.75, 0.875, 1)
poss.confounding.strength <- c("none", "mild", "strong")



for (ii in 1:nsims) {
  temp <- readRDS(paste0("Simulations_DAW/Replicate", ii, ".RDS"))
  
  for (aa in 1:length(poss.trt.prop)) {
    for (bb in 1:length(poss.trial.size)) {
      for (cc in 1:length(poss.ehr.fraction)) {
        for (dd in 1:length(poss.trt.effect)) {
          for (ff in 1:length(poss.confounding.strength)) {
            temp2 <- as.data.frame(temp[[aa]][[bb]][[cc]][[dd]][[ff]]$results)
            
            temp2$lb <- temp2$est-1.96*sqrt(temp2$var)
            temp2$ub <- temp2$est+1.96*sqrt(temp2$var)
            
            # For accessing correct true marginal values
            if (ff == 1) {gg <- 1}
            if (ff == 2) {gg <- 2}
            if (ff == 3) {gg <- 4}

            my.log.hrs[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- temp2$est
            my.hrs[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- exp(temp2$est)
            
            my.vars[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- temp2$var
            
            my.coverage.one.rounded[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- obtain_ci_coverage(round(temp2$lb,3), round(temp2$ub,3), log(1))
            my.coverage.truth.rounded[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- obtain_ci_coverage(round(temp2$lb,3), round(temp2$ub,3), round(marginal_HR[[dd]][[gg]],3))
            
            my.coverage.one[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- obtain_ci_coverage(temp2$lb, temp2$ub, log(1))
            my.coverage.truth[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- obtain_ci_coverage(temp2$lb, temp2$ub, marginal_HR[[dd]][[gg]])
            
            my.ESSs[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- temp2$ESS
            
            my.alphas[[aa]][[bb]][[cc]][[dd]][[ff]][ii, ] <- temp2$alpha
            
          }
        }
      }
    }
  }
  print(ii)
}



saveRDS(my.log.hrs, "Log_HRs.RDS")
saveRDS(my.hrs, "HRs.RDS")
saveRDS(my.vars, "Vars.RDS")
saveRDS(my.ESSs, "ESSs.RDS")
saveRDS(my.coverage.one.rounded, "Coverage_One_Rounded.RDS")
saveRDS(my.coverage.truth.rounded, "Coverage_Truth_Rounded.RDS")
saveRDS(my.coverage.one, "Coverage_One.RDS")
saveRDS(my.coverage.truth, "Coverage_Truth.RDS")
saveRDS(my.alphas, "Alphas.RDS")


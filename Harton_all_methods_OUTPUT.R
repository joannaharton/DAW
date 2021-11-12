# Joanna Harton
# Last updated: August 2021

library(plyr)


my.alphas <- readRDS("Alphas.RDS")
my.ests <- readRDS("Log_HRs.RDS")
my.vars <- readRDS("Vars.RDS")
my.ESSs <- readRDS("ESSs.RDS")
my.coverage.one <- readRDS("Coverage_One.RDS")
my.coverage.truth <- readRDS("Coverage_Truth.RDS")

# Estimated true marginal HR by running trial only method with a sample size of 1,000,000
true.marginal <- readRDS("true_marginal_HRs.RDS")

####### Plotting functions ########
my_bias_plot <- function(metric=c("actual", "percent"), trt.prop, trial.size, ehr.frac=1, confounding, my.y.min, my.y.max) {
  
  dat <- my.ests
  
  if (trt.prop==0.67) {aa <- 1}
  if (trt.prop==0.75) {aa <- 2}
  
  if (trial.size==100) {bb <- 1}
  if (trial.size==1000) {bb <- 2}
  
  if (ehr.frac==1) {cc <- 1}
  if (ehr.frac==5) {cc <- 2}
  
  # ff is for accessing correct true marginal values
  if (confounding=="none") {dd <- 1; ff <- 1}
  if (confounding=="mild") {dd <- 2; ff <- 2}
  if (confounding=="strong") {dd <- 3; ff <- 4}
  
  true.marg_0.5 <- exp(true.marginal$trt.0.5[[ff]])
  true.marg_0.75 <- exp(true.marginal$trt.0.75[[ff]])
  true.marg_0.875 <- exp(true.marginal$trt.0.875[[ff]])
  true.marg_1 <- exp(true.marginal$trt.1[[ff]])
  
  dat.temp <- dat[[aa]][[bb]][[cc]]
  
  dat_0.5 <- exp(apply(dat.temp[[1]][[dd]], MARGIN = 2, FUN = mean, na.rm=T))
  dat_0.75 <- exp(apply(dat.temp[[2]][[dd]], MARGIN = 2, FUN = mean, na.rm=T))
  dat_0.875 <- exp(apply(dat.temp[[3]][[dd]], MARGIN = 2, FUN = mean, na.rm=T))
  dat_1 <- exp(apply(dat.temp[[4]][[dd]], MARGIN = 2, FUN = mean, na.rm=T))
  
  if (metric=="percent") {
    dat_0.5 <- ((dat_0.5-true.marg_0.5)/true.marg_0.5)*100
    dat_0.75 <- ((dat_0.75-true.marg_0.75)/true.marg_0.75)*100
    dat_0.875 <- ((dat_0.875-true.marg_0.875)/true.marg_0.875)*100
    dat_1 <- ((dat_1-true.marg_1)/true.marg_1)*100
    
    new.dat <- rbind(dat_0.5, dat_0.75, dat_0.875, dat_1)
    
  } else {
    dat_0.5 <- dat_0.5-true.marg_0.5
    dat_0.75 <- dat_0.75-true.marg_0.75
    dat_0.875 <- dat_0.875-true.marg_0.875
    dat_1 <- dat_1-true.marg_1
    
    new.dat <- rbind(dat_0.5, dat_0.75, dat_0.875, dat_1)
    }
  
  new.dat <- as.data.frame(new.dat)
  
  

  my.trt.hrs <- c(true.marg_0.5, true.marg_0.75, true.marg_0.875, true.marg_1)
  my.main <- paste0("Trial Size: ", trial.size, "; Confounding Level: ", confounding)
  
  if (metric=="actual") {
    my.ylab <- "Bias"
  } else {my.ylab <- "Percent Bias"}
  
  
  # my.y.min <- min(new.dat)
  # my.y.max <- max(new.dat)
  
  my.x <- my.trt.hrs
  
  par(las=1)
  
  plot(my.x, new.dat$ie, type="b", pch=16, col="darkgreen", 
       ylim=c(my.y.min,my.y.max), xlim=c(min(my.x), max(my.x)), xaxt="n", ylab=my.ylab, xlab="Treatment HRs (Marginal)", main="")
  axis(1, at=my.x, labels=round(my.x, 3))
  abline(h=0, lty=3, lwd=2)
  points(my.x, new.dat$fp, type="b", pch=16, col="darkred")
  points(my.x, new.dat$pp_0.25, type="b", pch=16, col="deepskyblue")
  points(my.x, new.dat$pp_0.5, type="b", pch=16, col="dodgerblue")
  points(my.x, new.dat$pp_0.75, type="b", pch=16, col="dodgerblue4")
  points(my.x, new.dat$npp, type="b", pch=16, col="purple")
  
  points(my.x, new.dat$lin, type="b", pch=16, col="hotpink")
  
  points(my.x, new.dat$daw2, type="b", pch=16, col="darkorange")

  
}

my_variance_plot <- function(trt.prop, trial.size, ehr.frac=1, confounding) {
  
  dat <- my.vars
  
  if (trt.prop==0.67) {aa <- 1}
  if (trt.prop==0.75) {aa <- 2}
  
  if (trial.size==100) {bb <- 1}
  if (trial.size==1000) {bb <- 2}
  
  if (ehr.frac==1) {cc <- 1}
  if (ehr.frac==5) {cc <- 2}
  
  # ff is for accessing correct true marginal values
  if (confounding=="none") {dd <- 1; ff <- 1}
  if (confounding=="mild") {dd <- 2; ff <- 2}
  if (confounding=="strong") {dd <- 3; ff <- 4}
  
  true.marg_0.5 <- exp(true.marginal$trt.0.5[[ff]])
  true.marg_0.75 <- exp(true.marginal$trt.0.75[[ff]])
  true.marg_0.875 <- exp(true.marginal$trt.0.875[[ff]])
  true.marg_1 <- exp(true.marginal$trt.1[[ff]])
  
  dat.temp <- dat[[aa]][[bb]][[cc]]
  
  dat_0.5 <- apply(dat.temp[[1]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_0.75 <- apply(dat.temp[[2]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_0.875 <- apply(dat.temp[[3]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_1 <- apply(dat.temp[[4]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  

  new.dat <- rbind(dat_0.5, dat_0.75, dat_0.875, dat_1)
  new.dat <- as.data.frame(new.dat)
  
  my.trt.hrs <- c(true.marg_0.5, true.marg_0.75, true.marg_0.875, true.marg_1)
  my.main <- paste0("Trial Size: ", trial.size, "; Confounding Level: ", confounding)
  
  my.y.min <- 0
  my.y.max <- max(new.dat)
  
  my.x <- my.trt.hrs
  
  par(las=1)
  
  plot(my.x, new.dat$ie, type="b", pch=16, col="darkgreen", 
       ylim=c(my.y.min,my.y.max), xlim=c(min(my.x), max(my.x)), xaxt="n", ylab="Variance", xlab="Treatment HRs (Marginal)", main="")
  axis(1, at=my.x, labels=round(my.x, 3))
  abline(h=0, lty=3, lwd=2)
  points(my.x, new.dat$fp, type="b", pch=16, col="darkred")
  points(my.x, new.dat$pp_0.25, type="b", pch=16, col="deepskyblue")
  points(my.x, new.dat$pp_0.5, type="b", pch=16, col="dodgerblue")
  points(my.x, new.dat$pp_0.75, type="b", pch=16, col="dodgerblue4")
  points(my.x, new.dat$npp, type="b", pch=16, col="purple")
  
  points(my.x, new.dat$lin, type="b", pch=16, col="hotpink")
  
  points(my.x, new.dat$daw2, type="b", pch=16, col="darkorange")

  
}

my_coverage_plot <- function(trt.prop, trial.size, ehr.frac=1, confounding, my.y.min) {
  
  dat <- my.coverage.truth
  
  if (trt.prop==0.67) {aa <- 1}
  if (trt.prop==0.75) {aa <- 2}
  
  if (trial.size==100) {bb <- 1}
  if (trial.size==1000) {bb <- 2}
  
  if (ehr.frac==1) {cc <- 1}
  if (ehr.frac==5) {cc <- 2}
  
  # ff is for accessing correct true marginal values
  if (confounding=="none") {dd <- 1; ff <- 1}
  if (confounding=="mild") {dd <- 2; ff <- 2}
  if (confounding=="strong") {dd <- 3; ff <- 4}
  
  true.marg_0.5 <- exp(true.marginal$trt.0.5[[ff]])
  true.marg_0.75 <- exp(true.marginal$trt.0.75[[ff]])
  true.marg_0.875 <- exp(true.marginal$trt.0.875[[ff]])
  true.marg_1 <- exp(true.marginal$trt.1[[ff]])
  
  dat.temp <- dat[[aa]][[bb]][[cc]]
  
  dat_0.5 <- apply(dat.temp[[1]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_0.75 <- apply(dat.temp[[2]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_0.875 <- apply(dat.temp[[3]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_1 <- apply(dat.temp[[4]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  
  
  new.dat <- rbind(dat_0.5, dat_0.75, dat_0.875, dat_1)
  new.dat <- as.data.frame(new.dat)
  
  
  my.trt.hrs <- c(true.marg_0.5, true.marg_0.75, true.marg_0.875, true.marg_1)
  my.main <- paste0("Trial Size: ", trial.size, "; Confounding Level: ", confounding)
  
  
  my.y.max <- 1
  
  my.x <- my.trt.hrs
  
  par(las=1)
  
  plot(my.x, new.dat$ie, type="b", pch=16, col="darkgreen", 
       ylim=c(my.y.min,my.y.max), xlim=c(min(my.x), max(my.x)), xaxt="n", ylab="Coverage", xlab="Treatment HRs (Marginal)", main="")
  axis(1, at=my.x, labels=round(my.x, 3))
  abline(h=0.95, lty=3, lwd=2)
  points(my.x, new.dat$fp, type="b", pch=16, col="darkred")
  points(my.x, new.dat$pp_0.25, type="b", pch=16, col="deepskyblue")
  points(my.x, new.dat$pp_0.5, type="b", pch=16, col="dodgerblue")
  points(my.x, new.dat$pp_0.75, type="b", pch=16, col="dodgerblue4")
  points(my.x, new.dat$npp, type="b", pch=16, col="purple")
  
  points(my.x, new.dat$lin, type="b", pch=16, col="hotpink")
  
  points(my.x, new.dat$daw2, type="b", pch=16, col="darkorange")

}

my_power_plot <- function(trt.prop, trial.size, ehr.frac=1, confounding) {
  
  dat <- my.coverage.one
  
  if (trt.prop==0.67) {aa <- 1}
  if (trt.prop==0.75) {aa <- 2}
  
  if (trial.size==100) {bb <- 1}
  if (trial.size==1000) {bb <- 2}
  
  if (ehr.frac==1) {cc <- 1}
  if (ehr.frac==5) {cc <- 2}
  
  # ff is for accessing correct true marginal values
  if (confounding=="none") {dd <- 1; ff <- 1}
  if (confounding=="mild") {dd <- 2; ff <- 2}
  if (confounding=="strong") {dd <- 3; ff <- 4}
  
  true.marg_0.5 <- exp(true.marginal$trt.0.5[[ff]])
  true.marg_0.75 <- exp(true.marginal$trt.0.75[[ff]])
  true.marg_0.875 <- exp(true.marginal$trt.0.875[[ff]])
  true.marg_1 <- exp(true.marginal$trt.1[[ff]])
  
  dat.temp <- dat[[aa]][[bb]][[cc]]
  
  dat_0.5 <- apply(dat.temp[[1]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_0.75 <- apply(dat.temp[[2]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_0.875 <- apply(dat.temp[[3]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  dat_1 <- apply(dat.temp[[4]][[dd]], MARGIN = 2, FUN = mean, na.rm=T)
  
  
  new.dat <- rbind(dat_0.5, dat_0.75, dat_0.875)
  new.dat <- as.data.frame(new.dat)
  new.dat <- 1-new.dat
  
  
  my.trt.hrs <- c(true.marg_0.5, true.marg_0.75, true.marg_0.875)
  my.main <- paste0("Trial Size: ", trial.size, "; Confounding Level: ", confounding)
  
  my.y.min <- 0
  my.y.max <- 1
  
  my.x <- my.trt.hrs
  
  par(las=1)
  
  plot(my.x, new.dat$ie, type="b", pch=16, col="darkgreen", 
       ylim=c(my.y.min,my.y.max), xlim=c(min(my.x), max(my.x)), xaxt="n", ylab="Power", xlab="Treatment HRs (Marginal)", main="")
  axis(1, at=my.x, labels=round(my.x, 3))
  abline(h=0.80, lty=3, lwd=2)
  points(my.x, new.dat$fp, type="b", pch=16, col="darkred")
  points(my.x, new.dat$pp_0.25, type="b", pch=16, col="deepskyblue")
  points(my.x, new.dat$pp_0.5, type="b", pch=16, col="dodgerblue")
  points(my.x, new.dat$pp_0.75, type="b", pch=16, col="dodgerblue4")
  points(my.x, new.dat$npp, type="b", pch=16, col="purple")
  
  points(my.x, new.dat$lin, type="b", pch=16, col="hotpink")
  
  points(my.x, new.dat$daw2, type="b", pch=16, col="darkorange")

}

my_ESS_table <- function() {
  
  dat <- my.ESSs
  
  dat_trt067 <- dat$trt.prop.0.67
  dat_trt075 <- dat$trt.prop.0.75
  
  dat_trt067_trial100 <- dat_trt067$trial.100$ehr.frac.1
  dat_trt075_trial100 <- dat_trt075$trial.100$ehr.frac.1
  dat_trt067_trial1000 <- dat_trt067$trial.1000$ehr.frac.1
  dat_trt075_trial1000 <- dat_trt075$trial.1000$ehr.frac.1
  
  
  dat_trt067_trial100_confnone <- dat_trt067_trial100$trt.1$no.confounding
  dat_trt075_trial100_confnone <- dat_trt075_trial100$trt.1$no.confounding
  dat_trt067_trial1000_confnone <- dat_trt067_trial1000$trt.1$no.confounding
  dat_trt075_trial1000_confnone <- dat_trt075_trial1000$trt.1$no.confounding

  dat_trt067_trial100_confmild <- dat_trt067_trial100$trt.1$mild.confounding
  dat_trt075_trial100_confmild <- dat_trt075_trial100$trt.1$mild.confounding
  dat_trt067_trial1000_confmild <- dat_trt067_trial1000$trt.1$mild.confounding
  dat_trt075_trial1000_confmild <- dat_trt075_trial1000$trt.1$mild.confounding

  dat_trt067_trial100_confstrong <- dat_trt067_trial100$trt.1$strong.confounding
  dat_trt075_trial100_confstrong <- dat_trt075_trial100$trt.1$strong.confounding
  dat_trt067_trial1000_confstrong <- dat_trt067_trial1000$trt.1$strong.confounding
  dat_trt075_trial1000_confstrong <- dat_trt075_trial1000$trt.1$strong.confounding

  
  
  
  trt067_trial100_confnone <- apply(dat_trt067_trial100_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial100_confnone <- apply(dat_trt075_trial100_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt067_trial1000_confnone <- apply(dat_trt067_trial1000_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial1000_confnone <- apply(dat_trt075_trial1000_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)

  trt067_trial100_confmild <- apply(dat_trt067_trial100_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial100_confmild <- apply(dat_trt075_trial100_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt067_trial1000_confmild <-apply(dat_trt067_trial1000_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial1000_confmild <- apply(dat_trt075_trial1000_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)

  trt067_trial100_confstrong <- apply(dat_trt067_trial100_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial100_confstrong <- apply(dat_trt075_trial100_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt067_trial1000_confstrong <- apply(dat_trt067_trial1000_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial1000_confstrong <- apply(dat_trt075_trial1000_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)

  

  conf_none <- cbind(trt067_trial100_confnone, trt075_trial100_confnone,
                     trt067_trial1000_confnone, trt075_trial1000_confnone)

  conf_mild <- cbind(trt067_trial100_confmild, trt075_trial100_confmild,
                     trt067_trial1000_confmild, trt075_trial1000_confmild)

  conf_strong <- cbind(trt067_trial100_confstrong, trt075_trial100_confstrong,
                     trt067_trial1000_confstrong, trt075_trial1000_confstrong)

  all_conf <- cbind(conf_none, conf_mild, conf_strong)
  
  all_conf <- as.data.frame(all_conf)
  all_conf <- round(all_conf, 2)
  
  all_conf <- all_conf[1:8,]
  return(all_conf)
}

my_type1error_table <- function() {
  
  dat <- my.coverage.one
  
  dat_trt067 <- dat$trt.prop.0.67
  dat_trt075 <- dat$trt.prop.0.75
  
  dat_trt067_trial100 <- dat_trt067$trial.100$ehr.frac.1
  dat_trt075_trial100 <- dat_trt075$trial.100$ehr.frac.1
  dat_trt067_trial1000 <- dat_trt067$trial.1000$ehr.frac.1
  dat_trt075_trial1000 <- dat_trt075$trial.1000$ehr.frac.1
  
  
  dat_trt067_trial100_confnone <- dat_trt067_trial100$trt.1$no.confounding
  dat_trt075_trial100_confnone <- dat_trt075_trial100$trt.1$no.confounding
  dat_trt067_trial1000_confnone <- dat_trt067_trial1000$trt.1$no.confounding
  dat_trt075_trial1000_confnone <- dat_trt075_trial1000$trt.1$no.confounding
  
  dat_trt067_trial100_confmild <- dat_trt067_trial100$trt.1$mild.confounding
  dat_trt075_trial100_confmild <- dat_trt075_trial100$trt.1$mild.confounding
  dat_trt067_trial1000_confmild <- dat_trt067_trial1000$trt.1$mild.confounding
  dat_trt075_trial1000_confmild <- dat_trt075_trial1000$trt.1$mild.confounding
  
  dat_trt067_trial100_confstrong <- dat_trt067_trial100$trt.1$strong.confounding
  dat_trt075_trial100_confstrong <- dat_trt075_trial100$trt.1$strong.confounding
  dat_trt067_trial1000_confstrong <- dat_trt067_trial1000$trt.1$strong.confounding
  dat_trt075_trial1000_confstrong <- dat_trt075_trial1000$trt.1$strong.confounding
  
  
  
  trt067_trial100_confnone <- 1-apply(dat_trt067_trial100_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial100_confnone <- 1-apply(dat_trt075_trial100_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt067_trial1000_confnone <- 1-apply(dat_trt067_trial1000_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial1000_confnone <- 1-apply(dat_trt075_trial1000_confnone, MARGIN = 2, FUN = mean, na.rm = TRUE)
  
  trt067_trial100_confmild <- 1-apply(dat_trt067_trial100_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial100_confmild <- 1-apply(dat_trt075_trial100_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt067_trial1000_confmild <-1-apply(dat_trt067_trial1000_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial1000_confmild <- 1-apply(dat_trt075_trial1000_confmild, MARGIN = 2, FUN = mean, na.rm = TRUE)
  
  trt067_trial100_confstrong <- 1-apply(dat_trt067_trial100_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial100_confstrong <- 1-apply(dat_trt075_trial100_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt067_trial1000_confstrong <- 1-apply(dat_trt067_trial1000_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  trt075_trial1000_confstrong <- 1-apply(dat_trt075_trial1000_confstrong, MARGIN = 2, FUN = mean, na.rm = TRUE)
  
  
  conf_none <- cbind(trt067_trial100_confnone, trt075_trial100_confnone,
                     trt067_trial1000_confnone, trt075_trial1000_confnone)
  
  conf_mild <- cbind(trt067_trial100_confmild, trt075_trial100_confmild,
                     trt067_trial1000_confmild, trt075_trial1000_confmild)
  
  conf_strong <- cbind(trt067_trial100_confstrong, trt075_trial100_confstrong,
                       trt067_trial1000_confstrong, trt075_trial1000_confstrong)
  
  all_conf <- cbind(conf_none, conf_mild, conf_strong)
  
  all_conf <- as.data.frame(all_conf)
  all_conf <- round(all_conf, 4)
  all_conf <- all_conf[1:8,]
  
  return(all_conf)
  
}

######

write.csv(my_ESS_table(), file="ESS_Table.csv")
write.csv(my_type1error_table(), file="TypeI_Error_Table.csv")



my.legend.text <- c("Trial Only", "Full Pooling", "PP, 0.25", "PP, 0.5", "PP, 0.75", "NPP", "Lin", "DAW")
my.legend.cols <- c("darkgreen", "darkred", "deepskyblue", "dodgerblue", "dodgerblue4", "purple", "hotpink", "darkorange")

###### Actual Bias Plots ######
pdf("Actual Bias.pdf", width=9, height=8)

my.actual.bias.y.min <- -0.25
my.actual.bias.y.max <- 0.025

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))
my_bias_plot(metric="actual", trt.prop=0.67, trial.size=100, confounding="mild", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)

text(x=0.35, y=-0.125, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)


my_bias_plot(metric="actual", trt.prop=0.67, trial.size=100, confounding="strong", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)


legend(0.2, 0.1, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 0.1, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 0.1, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 0.1, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

 
polygon(x=c(0.1,0.1,0.9,0.9), y=c(0.075,0.125,0.125,0.075), xpd=NA)

my_bias_plot(metric="actual", trt.prop=0.67, trial.size=1000, confounding="mild", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)

text(x=0.35, y=-0.125, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=-0.4, "Mild Confounding", xpd=NA, cex=1.2)

my_bias_plot(metric="actual", trt.prop=0.67, trial.size=1000, confounding="strong", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)

text(x=0.8, y=-0.4, "Strong Confounding", xpd=NA, cex=1.2)


####


my.actual.bias.y.min <- -0.25
my.actual.bias.y.max <- 0.025

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))
my_bias_plot(metric="actual", trt.prop=0.75, trial.size=100, confounding="mild", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)

text(x=0.35, y=-0.125, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)


my_bias_plot(metric="actual", trt.prop=0.75, trial.size=100, confounding="strong", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)


legend(0.2, 0.1, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 0.1, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 0.1, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 0.1, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(0.075,0.125,0.125,0.075), xpd=NA)

my_bias_plot(metric="actual", trt.prop=0.75, trial.size=1000, confounding="mild", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)

text(x=0.35, y=-0.125, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=-0.4, "Mild Confounding", xpd=NA, cex=1.2)

my_bias_plot(metric="actual", trt.prop=0.75, trial.size=1000, confounding="strong", 
             my.y.min=my.actual.bias.y.min, my.y.max=my.actual.bias.y.max)

text(x=0.8, y=-0.4, "Strong Confounding", xpd=NA, cex=1.2)


dev.off()



###### Variance Plots ######
pdf("Variance.pdf", width=9, height=8)

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_variance_plot(trt.prop=0.67, trial.size=100, confounding="mild")
text(x=0.35, y=0.025, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_variance_plot(trt.prop=0.67, trial.size=100, confounding="strong")

legend(0.2, 0.07, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 0.07, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 0.07, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 0.07, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(0.075,0.065,0.065,0.075), xpd=NA)

my_variance_plot(trt.prop=0.67, trial.size=1000, confounding="mild")

text(x=0.35, y=0.0025, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=-0.0025, "Mild Confounding", xpd=NA, cex=1.2)

my_variance_plot(trt.prop=0.67, trial.size=1000, confounding="strong")


text(x=0.8, y=-0.0025, "Strong Confounding", xpd=NA, cex=1.2)


####


par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_variance_plot(trt.prop=0.75, trial.size=100, confounding="mild")
text(x=0.35, y=0.03, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_variance_plot(trt.prop=0.75, trial.size=100, confounding="strong")

legend(0.2, 0.08, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 0.08, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 0.08, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 0.08, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(0.075,0.085,0.085,0.075), xpd=NA)

my_variance_plot(trt.prop=0.75, trial.size=1000, confounding="mild")

text(x=0.35, y=0.003, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=-0.003, "Mild Confounding", xpd=NA, cex=1.2)

my_variance_plot(trt.prop=0.75, trial.size=1000, confounding="strong")


text(x=0.8, y=-0.003, "Strong Confounding", xpd=NA, cex=1.2)





dev.off()





###### Coverage Plots ######
pdf("Coverage.pdf", width=9, height=8)

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_coverage_plot(trt.prop=0.67, trial.size=100, confounding="mild", my.y.min = 0.55)

text(x=0.35, y=0.8, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_coverage_plot(trt.prop=0.67, trial.size=100, confounding="strong", my.y.min = 0.55)

legend(0.2, 1.1, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 1.1, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 1.1, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 1.1, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(1.05,1.15,1.15,1.05), xpd=NA)

my_coverage_plot(trt.prop=0.67, trial.size=1000, confounding="mild", my.y.min = 0)

text(x=0.35, y=0.5, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=-0.5, "Mild Confounding", xpd=NA, cex=1.2)


my_coverage_plot(trt.prop=0.67, trial.size=1000, confounding="strong", my.y.min = 0)

text(x=0.8, y=-0.5, "Strong Confounding", xpd=NA, cex=1.2)

####

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_coverage_plot(trt.prop=0.75, trial.size=100, confounding="mild", my.y.min = 0.55)

text(x=0.35, y=0.8, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_coverage_plot(trt.prop=0.75, trial.size=100, confounding="strong", my.y.min = 0.55)

legend(0.2, 1.1, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 1.1, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 1.1, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 1.1, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(1.05,1.15,1.15,1.05), xpd=NA)

my_coverage_plot(trt.prop=0.75, trial.size=1000, confounding="mild", my.y.min = 0)

text(x=0.35, y=0.5, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=-0.5, "Mild Confounding", xpd=NA, cex=1.2)


my_coverage_plot(trt.prop=0.75, trial.size=1000, confounding="strong", my.y.min = 0)

text(x=0.8, y=-0.5, "Strong Confounding", xpd=NA, cex=1.2)



dev.off()


###### ZOOMED IN Coverage Plots ######

pdf("Coverage ZOOMED IN.pdf", width=9, height=8)

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_coverage_plot(trt.prop=0.67, trial.size=100, confounding="mild", my.y.min = 0.90)

text(x=0.35, y=0.95, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_coverage_plot(trt.prop=0.67, trial.size=100, confounding="strong", my.y.min = 0.9)

legend(0.2, 1.025, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 1.025, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 1.025, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 1.025, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(1.0125,1.0375,1.0375,1.0125), xpd=NA)

my_coverage_plot(trt.prop=0.67, trial.size=1000, confounding="mild", my.y.min = 0.9)

text(x=0.35, y=0.95, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=0.85, "Mild Confounding", xpd=NA, cex=1.2)


my_coverage_plot(trt.prop=0.67, trial.size=1000, confounding="strong", my.y.min = 0.9)

text(x=0.8, y=0.85, "Strong Confounding", xpd=NA, cex=1.2)

####

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_coverage_plot(trt.prop=0.75, trial.size=100, confounding="mild", my.y.min = 0.90)

text(x=0.35, y=0.95, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_coverage_plot(trt.prop=0.75, trial.size=100, confounding="strong", my.y.min = 0.9)

legend(0.2, 1.025, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.4, 1.025, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.6, 1.025, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.8, 1.025, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.1,0.1,0.9,0.9), y=c(1.0125,1.0375,1.0375,1.0125), xpd=NA)

my_coverage_plot(trt.prop=0.75, trial.size=1000, confounding="mild", my.y.min = 0.9)

text(x=0.35, y=0.95, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.75, y=0.85, "Mild Confounding", xpd=NA, cex=1.2)


my_coverage_plot(trt.prop=0.75, trial.size=1000, confounding="strong", my.y.min = 0.9)

text(x=0.8, y=0.85, "Strong Confounding", xpd=NA, cex=1.2)


dev.off()

###### Power Plots ######
pdf("Power.pdf", width=9, height=8)

par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_power_plot(trt.prop=0.67, trial.size=100, confounding="mild")
text(x=0.4, y=0.5, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_power_plot(trt.prop=0.67, trial.size=100, confounding="strong")

legend(0.275, 1.25, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.425, 1.25, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.575, 1.25, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.725, 1.25, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.2,0.2,0.8,0.8), y=c(1.125,1.375,1.375,1.125), xpd=NA)

my_power_plot(trt.prop=0.67, trial.size=1000, confounding="mild")

text(x=0.4, y=0.5, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.7, y=-0.55, "Mild Confounding", xpd=NA, cex=1.2)

my_power_plot(trt.prop=0.67, trial.size=1000, confounding="strong")

text(x=0.75, y=-0.55, "Strong Confounding", xpd=NA, cex=1.2)

####


par(las=1, mfrow=c(2,2), mar=c(5,5,1,1), oma=c(4, 4, 6, 1))

my_power_plot(trt.prop=0.75, trial.size=100, confounding="mild")
text(x=0.4, y=0.5, "Trial Size: 100", srt=90, xpd=NA, cex=1.2)

my_power_plot(trt.prop=0.75, trial.size=100, confounding="strong")

legend(0.275, 1.25, legend=my.legend.text[c(1,2)],
       col=my.legend.cols[c(1,2)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.425, 1.25, legend=my.legend.text[c(3,4)],
       col=my.legend.cols[c(3,4)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.575, 1.25, legend=my.legend.text[c(5,6)],
       col=my.legend.cols[c(5,6)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")

legend(0.725, 1.25, legend=my.legend.text[c(7,8)],
       col=my.legend.cols[c(7,8)],
       lty=1, pch=16, xpd=NA, horiz=F, xjust = 0.5, yjust=0.5, bty="n")


polygon(x=c(0.2,0.2,0.8,0.8), y=c(1.125,1.375,1.375,1.125), xpd=NA)

my_power_plot(trt.prop=0.75, trial.size=1000, confounding="mild")

text(x=0.4, y=0.5, "Trial Size: 1,000", srt=90, xpd=NA, cex=1.2)
text(x=0.7, y=-0.55, "Mild Confounding", xpd=NA, cex=1.2)

my_power_plot(trt.prop=0.75, trial.size=1000, confounding="strong")

text(x=0.75, y=-0.55, "Strong Confounding", xpd=NA, cex=1.2)
dev.off()

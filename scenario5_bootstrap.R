### Set working directory
# setwd("...");

### Install required packages if necessary
#install.packages("BCEE");

### Load required packages;
require(BCEE);

### Change working directory for saving simulation results
#setwd("...");

p = 5;
trueATE = trueB = 1;
rho = 0.6;
var.list = c(paste("U",1:p,sep=""));
nrep = 1000;

coefs_GBCEE_TMLEv1 = numeric(nrep);
sd_GBCEE_TMLEv1 = numeric(nrep);
CP_GBCEE_TMLEv1 = numeric(nrep);

n = 1000;
{
 set.seed(61947919);
 for(i in 1:nrep)
 {
  ### simulate data
  U = rmvnorm(n = n, mean = rep(1, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
  U2 = t(apply(U, 1, function(x){c(x%*%t(x))}));
  X = rbinom(n, size = 1, p = plogis(-5 + U[,3] + U[,4] + U[,5] + rowSums(U2)));
  Y = rnorm(n, mean = trueB*X + 0.5*U[,1] + 0.5*U[,2] + 0.1*U[,3], sd = 1);

  ### BCEE v1 (omega = 500*n**0.5)
  results = GBCEE(X = X, Y = Y, U = U, omega = 1000*n**0.25, family.X = "binomial",
          niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1,
          var.comp = "bootstrap");
  coefs_GBCEE_TMLEv1[i] = results$beta[1];
  sd_GBCEE_TMLEv1[i] = results$stderr[1];
  CP_GBCEE_TMLEv1[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];

  print(i/nrep*100);
 }
 
 ### Record data
 file.name = paste("results_scenario5bootstrap_n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
 save(coefs_GBCEE_TMLEv1, sd_GBCEE_TMLEv1, CP_GBCEE_TMLEv1, file = file.name);
}

mean(coefs_GBCEE_TMLEv1);
sd(coefs_GBCEE_TMLEv1);
mean(sd_GBCEE_TMLEv1);
mean(CP_GBCEE_TMLEv1);















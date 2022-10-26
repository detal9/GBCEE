### Set working directory
# setwd("...");

### Install required packages if necessary
#install.packages("BCEE");

### Load required packages;
require(BCEE);

### Change working directory for saving simulation results
# setwd("...");

### Initializing objects

p = 5;
trueATE = trueB = 1;
rho = 0.6;
var.list = c(paste("U",1:p,sep=""));
nrep = 1000;
coefs_GBCEE100 = numeric(nrep);
coefs_GBCEEmodal = numeric(nrep);
coefs_GBCEE1000 = numeric(nrep);

sd_GBCEE100 = numeric(nrep);
sd_GBCEEmodal = numeric(nrep);
sd_GBCEE1000 = numeric(nrep);

included_GBCEE100 = matrix(, nrow = nrep, ncol = p);
included_GBCEEmodal = matrix(, nrow = nrep, ncol = p);
included_GBCEE1000 = matrix(, nrow = nrep, ncol = p);

CP_GBCEE100 = numeric(nrep);
CP_GBCEEmodal = numeric(nrep);
CP_GBCEE1000 = numeric(nrep);


pb = txtProgressBar(min = 1, max = nrep, initial = 1, style = 3, file = "", char = "*");
for(n in c(200, 1000))
{
 set.seed(61947919);
 print(n);
 for(i in 1:nrep)
 {
  print(i/nrep*100);
  setTxtProgressBar(pb, i);

  ### simulate data
  U = rmvnorm(n = n, mean = rep(1, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
  #U2 = matrix(apply(U, 2, function(x,y){x*y}, U), ncol = 25, byrow = FALSE);
  U2 = t(apply(U, 1, function(x){c(x%*%t(x))}));
  X = rbinom(n, size = 1, p = plogis(-5 + U[,3] + U[,4] + U[,5] + rowSums(U2)));
  Y = rnorm(n, mean = trueB*X + 0.5*U[,1] + 0.5*U[,2] + 0.1*U[,3], sd = 1);

  ### BCEE 100
  results = GBCEE(X = X, Y = Y, U = U, omega = 100*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEE100[i] = results$beta[1];
  sd_GBCEE100[i] = results$stderr[1];
  included_GBCEE100[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  CP_GBCEE100[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];

  ### BCEE 1000
  results = GBCEE(X = X, Y = Y, U = U, omega = 1000*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEE1000[i] = results$beta[1];
  sd_GBCEE1000[i] = results$stderr[1];
  included_GBCEE1000[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  CP_GBCEE1000[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];

  ### BCEE modal
  results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEEmodal[i] = results$models.Y[1,p+1];
  sd_GBCEEmodal[i] = results$models.Y[1,p+2];
  included_GBCEEmodal[i, ] = results$models.Y[1,1:p];
  CP_GBCEEmodal[i] = trueATE > coefs_GBCEEmodal[i] - 1.96*sd_GBCEEmodal[i] & trueATE < coefs_GBCEEmodal[i] + 1.96*sd_GBCEEmodal[i];

  print(data.frame(i, Sys.time()));
 }
 
 ### Record data
 file.name = paste("results_scenario5Supplemental_n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
 save(coefs_GBCEE100, sd_GBCEE100, included_GBCEE100, CP_GBCEE100,
      coefs_GBCEEmodal, sd_GBCEEmodal, included_GBCEEmodal, CP_GBCEEmodal,
      coefs_GBCEE1000, sd_GBCEE1000, included_GBCEE1000, CP_GBCEE1000,
      file = file.name);
}

















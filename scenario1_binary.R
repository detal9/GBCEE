### Set working directory
# setwd("...");

### Install required packages if necessary
#install.packages("mvtnorm");
#install.packages("BMA");
#install.packages("rms");
#install.packages("ctmle");
#install.packages("bacr");
#install.packages("madr");
#install.packages("lars"); # Required for BayesPen
#install.packages("SuppDists"); # Required for BayesPen
#install.packages("hdm");
#install.packages("penalized");
# BayesPen needs to be downloaded here: https://cran.r-project.org/src/contrib/Archive/BayesPen/
#  and then manually installed

### Load required packages;
require(MASS);
require(mvtnorm);
require(BMA);
require(ctmle);
require(boot);
require(bacr);
require(BCEE);
require(penalized);
source("OAL_v3.R");



### Change working directory for saving simulation results
setwd("...");

### Initializing objects

p = 40;
trueB = 1;
var.list = c(paste("U",1:p,sep=""));
nrep = 1000;
#n = 1000000;
#U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1, nrow = p) + matrix(0, nrow = p, ncol = p));
#U[,1:5] = 1*rowSums(U[,11:15]) + matrix(rnorm(n*5), ncol = 5);
#Y0 = rbinom(n, size = 1, p = plogis(0.1*rowSums(U[,1:10]) + trueB*0));
#Y1 = rbinom(n, size = 1, p = plogis(0.1*rowSums(U[,1:10]) + trueB*1));
#trueATE = mean(Y1) - mean(Y0);
trueATE = 0.189351;


coefs_GBCEE_TMLEv2 = numeric(nrep);
coefs_CTMLE = numeric(nrep);
coefs_BAC = numeric(nrep);
coefs_crude = numeric(nrep);
coefs_full = numeric(nrep);
coefs_target = numeric(nrep);
coefs_full2 = numeric(nrep);
coefs_target2 = numeric(nrep);

sd_GBCEE_TMLEv2 = numeric(nrep);
sd_CTMLE = numeric(nrep);
sd_BAC = numeric(nrep);
sd_crude = numeric(nrep);
sd_full = numeric(nrep);
sd_target = numeric(nrep);
sd_full2 = numeric(nrep);
sd_target2 = numeric(nrep);

included_GBCEEv2 = matrix(, nrow = nrep, ncol = p);
included_GBCEEXv2 = matrix(, nrow = nrep, ncol = p);
included_CTMLE = matrix(, nrow = nrep, ncol = p);
included_BAC = matrix(, nrow = nrep, ncol = p);



CP_GBCEE_TMLEv2 = numeric(nrep);
CP_CTMLE = numeric(nrep);
CP_GLIDER = numeric(nrep);
CP_BAC = numeric(nrep);
CP_crude = numeric(nrep);
CP_full = numeric(nrep);
CP_target = numeric(nrep);
CP_full2 = numeric(nrep);
CP_target2 = numeric(nrep);

coefs_OAL = sd_OAL = CP_OAL = numeric(nrep);
coefs_OAL_TMLE = sd_OAL_TMLE = CP_OAL_TMLE = numeric(nrep)
included_OAL = matrix(, nrow = nrep, ncol = p);


pb = txtProgressBar(min = 1, max = nrep, initial = 1, style = 3, file = "", char = "*");
for(n in c(1000))
{
 set.seed(61947919);
 for(i in 1:nrep)
 {
  setTxtProgressBar(pb, i);

  ### simulate data
  U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1, nrow = p) + matrix(0, nrow = p, ncol = p));
  U[,1:5] = 1*rowSums(U[,11:15]) + matrix(rnorm(n*5), ncol = 5);
  X = rbinom(n, size = 1, p = plogis(1*rowSums(U[,c(11:30)])));
  Y = rbinom(n, size = 1, p = plogis(0.1*rowSums(U[,1:10]) + trueB*X));


  ### BCEE v2 (omega = 500*n**0.5)
  results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "binomial", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEE_TMLEv2[i] = results$beta[1];
  sd_GBCEE_TMLEv2[i] = results$stderr[1];
  included_GBCEEv2[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  included_GBCEEXv2[i, ] = colSums(results$models.X[,1:p]*results$models.X[,p+1]);
  CP_GBCEE_TMLEv2[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];

  ### Outcome-adaptive lasso
  results = OAL.logistic(Y = Y, A = X, X = U);
  coefs_OAL[i] = results$ATE;
  coefs_OAL_TMLE[i] = results$TMLE;
  included_OAL[i,] = results$included;
 
  ### C-TMLE
  results = ctmleDiscrete(Y = Y, A = X, W = data.frame(U), preOrder = FALSE, family = "binomial");
  coefs_CTMLE[i] = results$est;
  sd_CTMLE[i] = sqrt(results$var.psi);
  included_CTMLE[i, ] = paste("X", seq(1:p), sep = "") %in% head(results$candidates$terms, results$best_k);
  CP_CTMLE[i] = trueATE >  coefs_CTMLE[i] - 1.96*sd_CTMLE[i] & trueATE <  coefs_CTMLE[i] + 1.96*sd_CTMLE[i];

  ### BAC
  ds = data.frame(Y,X,U);
  colnames(ds) = c("Y", "X", paste("U", 1:p, sep = ""));
  varlist = paste("U", 1:p, sep = "");
  results = bac(data = ds, exposure = "X", outcome = "Y", confounders = varlist, interactors = NULL,
                familyX = "binomial", family = "binomial", num_its = 5000, burnM = 500, burnB = 500, thin = 5);
  coefs_BAC[i] = mean(results$ACE);
  sd_BAC[i] = sd(results$ACE);
  included_BAC[i,] = summary(results)$PIP;
  CP_BAC[i] = trueATE < quantile(results$ACE, 0.975) & trueATE > quantile(results$ACE, 0.025);


  ### Crude
  crude = glm(Y~X, family = binomial(link = "identity"));
  coefs_crude[i] = coef(crude)[2];
  sd_crude[i] = sqrt(vcov(crude)[2,2]);
  CP_crude[i] = trueATE > coefs_crude[i] - 1.96*sd_crude[i] & trueATE < coefs_crude[i] + 1.96*sd_crude[i];

  ### Fully adjusted model
  full = glm(Y~X+U, family = "binomial");
  coefs_full[i] = mean(predict(full, newdata = data.frame(X = 1, ds[,-c(1,2)]), type = "res") - 
                predict(full, newdata = data.frame(X = 0, ds[,-c(1,2)]), type = "res"));


  ### Target model
  target = glm(Y~X+U[,1:10], family = "binomial");
  coefs_target[i] = mean(predict(target, newdata = data.frame(X = 1, ds[,-c(1,2)]), type = "res") - 
                predict(target, newdata = data.frame(X = 0, ds[,-c(1,2)]), type = "res"));

  ### Fully adjusted model - AIPTW
  fullX = glm(X~U, family = "binomial");
  pX1 = plogis(cbind(1, U)%*%coef(fullX));
  bounds = quantile(pX1, c(0.01, 0.99)); 
  pX1 = pmin(pX1, bounds[2]);
  pX1 = pmax(pX1, bounds[1]);
  H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
  Q1_0 = plogis(cbind(1, 1, U)%*%coef(full));
  Q0_0 = plogis(cbind(1, 0, U)%*%coef(full));
  QX_0 = plogis(cbind(1, X, U)%*%coef(full));
  coefs_full2[i] = mean((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0);
  sd_full2[i] = sqrt(var((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0 - coefs_full2[i])/n);
  CP_full2[i] = trueATE > coefs_full2[i] - 1.96*sd_full2[i] & trueATE < coefs_full2[i] + 1.96*sd_full2[i];

  ### Target model - AIPTW
  targetX = glm(X~U[,1:10], family = "binomial");
  pX1 = plogis(cbind(1, U[,1:10])%*%coef(targetX));
  bounds = quantile(pX1, c(0.01, 0.99)); 
  pX1 = pmin(pX1, bounds[2]);
  pX1 = pmax(pX1, bounds[1]);
  H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
  Q1_0 = plogis(cbind(1, 1, U[,1:10])%*%coef(target));
  Q0_0 = plogis(cbind(1, 0, U[,1:10])%*%coef(target));
  QX_0 = plogis(cbind(1, X, U[,1:10])%*%coef(target));
  coefs_target2[i] = mean((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0);
  sd_target2[i] = sqrt(var((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0 - coefs_target2[i])/n);
  CP_target2[i] = trueATE > coefs_target2[i] - 1.96*sd_target2[i] & trueATE < coefs_target2[i] + 1.96*sd_target2[i];

  print(data.frame(i, Sys.time()));
 }
 
 ### Record data
 file.name = paste("results_scenario1B_n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
 save(coefs_GBCEE_TMLEv2, 
      coefs_CTMLE, coefs_BAC, 
      coefs_crude, coefs_full, coefs_target, coefs_full2, coefs_target2,
      sd_GBCEE_TMLEv2, 
      sd_CTMLE, sd_BAC, sd_crude,
      sd_full, sd_target, sd_full2, sd_target2,
      included_GBCEEv2, included_GBCEEXv2,
      included_CTMLE, included_BAC, 
      CP_GBCEE_TMLEv2,
      CP_CTMLE, CP_BAC, CP_crude,
      CP_full, CP_target, CP_full2, CP_target2,
      coefs_OAL, coefs_OAL_TMLE, sd_OAL, sd_OAL_TMLE,
      included_OAL, CP_OAL, CP_OAL_TMLE,
      file = file.name);
}















### Set working directory
# setwd("...");

### Install required packages if necessary
#install.packages("penalized"); 
#install.packages("mvtnorm");
#install.packages("BMA");
#install.packages("rms");
#install.packages("ctmle");
#install.packages("bacr");
#install.packages("madr");
#install.packages("lars"); # Required for BayesPen
#install.packages("SuppDists"); # Required for BayesPen
#install.packages("hdm");
#install.packages("tmle");
# BayesPen needs to be downloaded here: https://cran.r-project.org/src/contrib/Archive/BayesPen/
#  and then manually installed

### Load required packages;
require(MASS);
require(penalized);
require(mvtnorm);
require(BMA);
require(ctmle);
require(boot);
require(bacr);
require(madr);
require(BayesPen);
require(hdm);
require(BCEE);
require(tmle);
source("OAL_V3.R");
source("glider.R");

library(devtools)
install_github(repo = "jantonelli111/HDconfounding")
library(HDconfounding)

### Change working directory for saving simulation results
# setwd("...");

### Initializing objects

p = 20;
rho = 0.5;
trueATE = trueB = 2;
beta = matrix(c(rep(0.6, 4), 0, 0,  rep(0, p - 6)), ncol = 1);
nu = matrix(c(1, 1, 0, 0, 1, 1, rep(0, p - 6)), ncol = 1);
var.list = c(paste("U",1:p,sep=""));
nrep = 1000;


coefs_GBCEE_TMLEv2 = numeric(nrep);
coefs_BAC = numeric(nrep);
coefs_BP = numeric(nrep);
coefs_crude = numeric(nrep);
coefs_full = numeric(nrep);
coefs_target = numeric(nrep);
coefs_full2 = numeric(nrep);
coefs_target2 = numeric(nrep);
coefs_SS = numeric(nrep);

sd_GBCEE_TMLEv2 = numeric(nrep);
sd_BAC = numeric(nrep);
sd_BP = numeric(nrep);
sd_crude = numeric(nrep);
sd_full = numeric(nrep);
sd_target = numeric(nrep);
sd_full2 = numeric(nrep);
sd_target2 = numeric(nrep);
sd_SS = numeric(nrep);

included_GBCEEv2 = matrix(, nrow = nrep, ncol = p);
included_GBCEEXv2 = matrix(, nrow = nrep, ncol = p);
included_BAC = matrix(, nrow = nrep, ncol = p);
included_BP = matrix(, nrow = nrep, ncol = p);
included_SS = matrix(, nrow = nrep, ncol = p);

CP_GBCEE_TMLEv2 = numeric(nrep);
CP_BAC = numeric(nrep);
CP_BP = numeric(nrep);
CP_crude = numeric(nrep);
CP_full = numeric(nrep);
CP_target = numeric(nrep);
CP_full2 = numeric(nrep);
CP_target2 = numeric(nrep);
CP_SS = numeric(nrep);

pb = txtProgressBar(min = 1, max = nrep, initial = 1, style = 3, file = "", char = "*");
for(n in c(200, 1000))
{
 set.seed(61947919);
 print(n);
 for(i in 1:nrep)
 {
  setTxtProgressBar(pb, i);

  ### simulate data
  U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
  X = rnorm(n, mean = U%*%nu);
  Y = rnorm(n, mean = trueB*X + U%*%beta, sd = 1);

  ### BCEE v2 (omega = 500*n**0.5)
  results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "gaussian",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEE_TMLEv2[i] = results$beta[1];
  sd_GBCEE_TMLEv2[i] = results$stderr[1];
  included_GBCEEv2[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  included_GBCEEXv2[i, ] = colSums(results$models.X[,1:p]*results$models.X[,p+1]);
  CP_GBCEE_TMLEv2[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];


  ### BAC
  ds = data.frame(Y,X,U);
  colnames(ds) = c("Y", "X", paste("U", 1:p, sep = ""));
  varlist = paste("U", 1:p, sep = "");
  results = bac(data = ds, exposure = "X", outcome = "Y", confounders = varlist, interactors = NULL,
                familyX = "gaussian", family = "gaussian", num_its = 5000, burnM = 500, burnB = 500, thin = 5);
  coefs_BAC[i] = mean(results$ACE);
  sd_BAC[i] = sd(results$ACE);
  included_BAC[i,] = summary(results)$PIP;
  CP_BAC[i] = trueATE < quantile(results$ACE, 0.975) & trueATE > quantile(results$ACE, 0.025);


  ### BayesPen
  fit1 = lm(Y ~ X + U);
  betahat = coef(fit1);
  cov = vcov(fit1);
  fit2 = glm(X ~ U, family = "gaussian");
  gammahat = coef(fit2);
  fit.BayesPen = BayesPen(beta = betahat, beta_cov = cov, confounder.weights = c(0, gammahat), force = c(1,2));
  path = fit.BayesPen$joint.path[-1,-c(1,2)];
  lag.path = rbind(fit.BayesPen$joint.path[-nrow(fit.BayesPen$joint.path),-c(1,2)]);
  to.add = (path - lag.path)>0;
  varnames = as.character(colnames(to.add));
  minalpha = 0;
  j = 1;
  varlist = NULL;
  ds = data.frame(Y,X,U);
  colnames(ds) = c("Y", "X", varnames);
  while(minalpha <= 0.25 & j <= nrow(to.add))
  {
   varlist.temp = c(varlist,varnames[to.add[j,] == 1]);
   y.form = formula(paste("Y~X+",paste(varlist.temp,collapse="+")));
   x.form = formula(paste("X~",paste(varlist.temp,collapse="+")));
   alpha1 = tail(summary(lm(y.form, data = ds))$coefficients[,4],1);
   alpha2 = tail(summary(glm(x.form, data = ds, family = "gaussian"))$coefficients[,4],1);
   minalpha = min(alpha1, alpha2);
   j = j + 1;
   if(minalpha <= 0.25)
   {
    varlist = varlist.temp;
   }
  }
  y.form = formula(paste("Y~X+",paste(varlist,collapse="+")));
  results = lm(y.form, data = ds);
  coefs_BP[i] = coef(results)[2];
  sd_BP[i] = sqrt(vcov(results)[2,2]);
  CP_BP[i] = trueATE >  coefs_BP[i] - 1.96*sd_BP[i] & trueATE <  coefs_BP[i] + 1.96*sd_BP[i];
  included_BP[i, ] = paste("U", seq(1:p), sep = "") %in% varlist;


  ### Crude
  crude = lm(Y~X);
  coefs_crude[i] = coef(crude)[2];
  sd_crude[i] = sqrt(vcov(crude)[2,2]);
  CP_crude[i] = trueATE > coefs_crude[i] - 1.96*sd_crude[i] & trueATE < coefs_crude[i] + 1.96*sd_crude[i];


  ### Fully adjusted model
  full = lm(Y~X+U);
  coefs_full[i] = coef(full)[2];
  sd_full[i] = sqrt(vcov(full)[2,2]);
  CP_full[i] = trueATE > coefs_full[i] - 1.96*sd_full[i] & trueATE < coefs_full[i] + 1.96*sd_full[i];


  ### Target model
  target = lm(Y~X+U[,1:4]);
  coefs_target[i] = coef(target)[2];
  sd_target[i] = sqrt(vcov(target)[2,2]);
  CP_target[i] = trueATE > coefs_target[i] - 1.96*sd_target[i] & trueATE < coefs_target[i] + 1.96*sd_target[i];


  ### Fully adjusted model - AIPTW
  fullX = lm(X~U);
  B0 = coef(full)[2]; # Initial estimate of slope
  Q0 = predict(full); # Y hat
  r0 = Q0 - B0*X; 
  predX0 = predict(fullX); # X hat
  Hg = X - predX0; # Clever covariate
  epsilon = coef(lm(Y~-1+Hg, offset = Q0)); # Fluctuation of estimate
  coefs_full2[i] = B1 = B0 + epsilon; # Final TMLE estimate of slope
  r1 = r0 - epsilon*predX0; 
  Q1 = B1*X + r1; # Updated Y hat
  k = -mean(Hg*X); 
  D = Hg*(Y - Q1); # Score
  IC = 1/k*D; # Influence curve
  sd_full2[i] = sqrt(mean(IC**2)/n);
  CP_full2[i] = trueATE > coefs_full2[i] - 1.96*sd_full2[i] & trueATE < coefs_full2[i] + 1.96*sd_full2[i];

  ### Target model - AIPTW
  targetX = lm(X~U[,1:4]);
  B0 = coef(target)[2]; # Initial estimate of slope
  Q0 = predict(target); # Y hat
  r0 = Q0 - B0*X; 
  predX0 = predict(targetX); # X hat
  Hg = X - predX0; # Clever covariate
  epsilon = coef(lm(Y~-1+Hg, offset = Q0)); # Fluctuation of estimate
  coefs_target2[i] = B1 = B0 + epsilon; # Final TMLE estimate of slope
  r1 = r0 - epsilon*predX0; 
  Q1 = B1*X + r1; # Updated Y hat
  k = -mean(Hg*X); 
  D = Hg*(Y - Q1); # Score
  IC = 1/k*D; # Influence curve
  sd_target2[i] = sqrt(mean(IC**2)/n);
  CP_target2[i] = trueATE > coefs_target2[i] - 1.96*sd_target2[i] & trueATE < coefs_target2[i] + 1.96*sd_target2[i];

  ## Spike and Slab
  tryCatch({
  sslEB = SSL(y=Y, z=X, x=U, nScans=20000, burn=10000, thin=10, lambda0="EB", z_type = "continuous");
  coefs_SS[i] = sslEB$TreatEffect;
  sd_SS[i] = sd(sslEB$TreatEffectPost);
  included_SS[i,] = sslEB$gammaPostMean;
  CP_SS[i] = trueATE > sslEB$TreatEffectCI[1] & trueATE < sslEB$TreatEffectCI[2];},
   error = function(e){
    coefs_SS[i] = NA; sd_SS[i] = NA; CP_SS[i] = NA;});

  print(i);
 }
 
 ### Record data
 file.name = paste("results_scenario2CC_n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
 save(coefs_GBCEE_TMLEv2, coefs_BAC, coefs_BP, coefs_crude, coefs_full, coefs_target, coefs_full2, coefs_target2,
      sd_GBCEE_TMLEv2, sd_BAC, sd_BP, sd_crude, sd_full, sd_target, sd_full2, sd_target2,
      CP_GBCEE_TMLEv2, CP_BAC, sd_BP, CP_crude, CP_full, CP_target, CP_full2, CP_target2,
      included_GBCEEv2, included_GBCEEXv2, included_BAC, included_BP,
      coefs_SS, sd_SS, included_SS, CP_SS, 
      file = file.name);
}
















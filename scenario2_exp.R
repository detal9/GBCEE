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
coefs_CTMLE = numeric(nrep);
coefs_GLIDER = numeric(nrep);
coefs_BAC = numeric(nrep);
coefs_MADR = numeric(nrep);
coefs_BP = numeric(nrep);
coefs_HDM = numeric(nrep);
coefs_crude = numeric(nrep);
coefs_full = numeric(nrep);
coefs_target = numeric(nrep);
coefs_full2 = numeric(nrep);
coefs_target2 = numeric(nrep);
coefs_SS = numeric(nrep);

sd_GBCEE_TMLEv2 = numeric(nrep);
sd_CTMLE = numeric(nrep);
sd_GLIDER = numeric(nrep);
sd_BAC = numeric(nrep);
sd_MADR = numeric(nrep);
sd_BP = numeric(nrep);
sd_HDM = numeric(nrep);
sd_crude = numeric(nrep);
sd_full = numeric(nrep);
sd_target = numeric(nrep);
sd_full2 = numeric(nrep);
sd_target2 = numeric(nrep);
sd_SS = numeric(nrep);


included_GBCEEv2 = matrix(, nrow = nrep, ncol = p);
included_GBCEEXv2 = matrix(, nrow = nrep, ncol = p);
included_CTMLE = matrix(, nrow = nrep, ncol = p);
included_GLIDER = matrix(, nrow = nrep, ncol = p);
included_BAC = matrix(, nrow = nrep, ncol = p);
included_MADRx = matrix(, nrow = nrep, ncol = p);
included_MADRy = matrix(, nrow = nrep, ncol = p);
included_BP = matrix(, nrow = nrep, ncol = p);
included_SS = matrix(, nrow = nrep, ncol = p);

CP_GBCEE_TMLEv2 = numeric(nrep);
CP_CTMLE = numeric(nrep);
CP_GLIDER = numeric(nrep);
CP_BAC = numeric(nrep);
CP_MADR = numeric(nrep);
CP_BP = numeric(nrep);
CP_HDM = numeric(nrep);
CP_crude = numeric(nrep);
CP_full = numeric(nrep);
CP_target = numeric(nrep);
CP_full2 = numeric(nrep);
CP_target2 = numeric(nrep);
CP_SS = numeric(nrep);

coefs_OAL = sd_OAL = CP_OAL = numeric(nrep);
coefs_OAL_TMLE = sd_OAL_TMLE = CP_OAL_TMLE = numeric(nrep)
included_OAL = matrix(, nrow = nrep, ncol = p);


stat.madr = function(x, i)
{
 Y = as.matrix(ds$Y[i]);
 X = as.matrix(ds$X[i]);
 U = as.matrix(ds[i,-c(1,2)]);
 madr(Y = Y, X = X, U = U)$madr;
}


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
  X = rbinom(n, 1, p = plogis(U%*%nu));
  Y = trueB*X + U%*%beta + rexp(n) - 1;

  ### BCEE v2 (omega = 500*n**0.5)
  tryCatch({
  results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEE_TMLEv2[i] = results$beta[1];
  sd_GBCEE_TMLEv2[i] = results$stderr[1];
  included_GBCEEv2[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  included_GBCEEXv2[i, ] = colSums(results$models.X[,1:p]*results$models.X[,p+1]);
  CP_GBCEE_TMLEv2[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];}, error = function(e){});

  ### Outcome-adaptive lasso
  tryCatch({
  results = OAL.lm(Y = Y, A = X, X = U);
  coefs_OAL[i] = results$ATE;
  coefs_OAL_TMLE[i] = results$TMLE;
  included_OAL[i,] = results$included;}, error = function(e){});
 
  ### C-TMLE
  tryCatch({
  results = ctmleDiscrete(Y = Y, A = X, W = data.frame(U), preOrder = FALSE);
  coefs_CTMLE[i] = results$est;
  sd_CTMLE[i] = sqrt(results$var.psi);
  included_CTMLE[i, ] = paste("X", seq(1:p), sep = "") %in% head(results$candidates$terms, results$best_k);
  CP_CTMLE[i] = trueATE >  coefs_CTMLE[i] - 1.96*sd_CTMLE[i] & trueATE <  coefs_CTMLE[i] + 1.96*sd_CTMLE[i];}, error = function(e){});

  ### GLIDER
  tryCatch({
  results = GLiDeR(Xorig = U, Yorig = Y, Aorig = X);
  coefs_GLIDER[i] = results$delta;
  included_GLIDER[i, ] = abs(results$alpha) > 0.00001;}, error = function(e){});

  ### BAC
  ds = data.frame(Y,X,U);
  colnames(ds) = c("Y", "X", paste("U", 1:p, sep = ""));
  varlist = paste("U", 1:p, sep = "");
  tryCatch({
  results = bac(data = ds, exposure = "X", outcome = "Y", confounders = varlist, interactors = NULL,
                familyX = "binomial", family = "gaussian", num_its = 5000, burnM = 500, burnB = 500, thin = 5);
  coefs_BAC[i] = mean(results$ACE);
  sd_BAC[i] = sd(results$ACE);
  included_BAC[i,] = summary(results)$PIP;
  CP_BAC[i] = trueATE < quantile(results$ACE, 0.975) & trueATE > quantile(results$ACE, 0.025);},
  error = function(e){});


  ### MADR
  tryCatch({
  results = madr(Y = Y, X = X, U = U);
  coefs_MADR[i] = results$madr;
  included_MADRx[i, ] = results$weight.ps;
  included_MADRy[i, ] = results$weight.om;}, error = function(e){});


  ### BayesPen
  tryCatch({
  fit1 = lm(Y ~ X + U);
  betahat = coef(fit1);
  cov = vcov(fit1);
  fit2 = glm(X ~ U, family = "binomial");
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
   alpha2 = tail(summary(glm(x.form, data = ds, family = "binomial"))$coefficients[,4],1);
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
  included_BP[i, ] = paste("U", seq(1:p), sep = "") %in% varlist;}, error = function(e){});


  ### HDM
  tryCatch({
  results = rlassoEffect(x = U, y = Y, d = X);
  coefs_HDM[i] = results$coefficients;
  sd_HDM[i] = results$se;
  CP_HDM[i] = trueATE >  coefs_HDM[i] - 1.96*sd_HDM[i] & trueATE <  coefs_HDM[i] + 1.96*sd_HDM[i];}, error = function(e){});


  ### Crude
  crude = lm(Y~X);
  coefs_crude[i] = coef(crude)[2];
  sd_crude[i] = sqrt(vcov(crude)[2,2]);
  CP_crude[i] = trueATE > coefs_crude[i] - 1.96*sd_crude[i] & trueATE < coefs_crude[i] + 1.96*sd_crude[i];

  ### Fully adjusted model
  tryCatch({
  full = lm(Y~X+U);
  coefs_full[i] = coef(full)[2];
  sd_full[i] = sqrt(vcov(full)[2,2]);
  CP_full[i] = trueATE > coefs_full[i] - 1.96*sd_full[i] & trueATE < coefs_full[i] + 1.96*sd_full[i];}, error = function(e){});


  ### Target model
  tryCatch({
  target = lm(Y~X+U[,1:4]);
  coefs_target[i] = coef(target)[2];
  sd_target[i] = sqrt(vcov(target)[2,2]);
  CP_target[i] = trueATE > coefs_target[i] - 1.96*sd_target[i] & trueATE < coefs_target[i] + 1.96*sd_target[i];}, error = function(e){});


  ### Fully adjusted model - AIPTW
  tryCatch({
  fullX = glm(X~U, family = "binomial");
  pX1 = plogis(cbind(1, U)%*%coef(fullX));
  bounds = quantile(pX1, c(0.01, 0.99)); 
  pX1 = pmin(pX1, bounds[2]);
  pX1 = pmax(pX1, bounds[1]);
  H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
  Q1_0 = cbind(1, 1, U)%*%coef(full);
  Q0_0 = cbind(1, 0, U)%*%coef(full);
  QX_0 = cbind(1, X, U)%*%coef(full);
  coefs_full2[i] = mean((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0);
  sd_full2[i] = sqrt(var((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0 - coefs_full2[i])/n);
  CP_full2[i] = trueATE > coefs_full2[i] - 1.96*sd_full2[i] & trueATE < coefs_full2[i] + 1.96*sd_full2[i];}, error = function(e){});


  ### Target model - AIPTW
  tryCatch({
  targetX = lm(X~U[,1:4]);
  pX1 = plogis(cbind(1, U[,1:4])%*%coef(targetX));
  bounds = quantile(pX1, c(0.01, 0.99)); 
  pX1 = pmin(pX1, bounds[2]);
  pX1 = pmax(pX1, bounds[1]);
  H1 = X/pX1; H0 = (1 - X)/(1 - pX1);
  Q1_0 = cbind(1, 1, U[,1:4])%*%coef(target);
  Q0_0 = cbind(1, 0, U[,1:4])%*%coef(target);
  QX_0 = cbind(1, X, U[,1:4])%*%coef(target);
  coefs_target2[i] = mean((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0);
  sd_target2[i] = sqrt(var((H1 - H0)*(Y - QX_0) + Q1_0 - Q0_0 - coefs_target2[i])/n);
  CP_target2[i] = trueATE > coefs_target2[i] - 1.96*sd_target2[i] & trueATE < coefs_target2[i] + 1.96*sd_target2[i];}, error = function(e){});

  sslEB = SSL(y=Y, z=X, x=U, nScans=20000, burn=10000, thin=10, lambda0="EB");
  coefs_SS[i] = sslEB$TreatEffect;
  sd_SS[i] = sd(sslEB$TreatEffectPost);
  included_SS[i,] = sslEB$gammaPostMean;
  CP_SS[i] = trueATE > sslEB$TreatEffectCI[1] & trueATE < sslEB$TreatEffectCI[2];
 }
 
 ### Record data
 file.name = paste("results_scenario2Exp_n",n,"_",gsub("-","",Sys.Date()),".Rdata", sep="");
 save(coefs_GBCEE_TMLEv2, 
      coefs_OAL, coefs_CTMLE, coefs_GLIDER, coefs_BAC, coefs_MADR, coefs_BP, coefs_HDM, 
      coefs_crude, coefs_full, coefs_target, coefs_full2, coefs_target2,
      sd_GBCEE_TMLEv2, 
      sd_OAL, sd_CTMLE, sd_GLIDER, sd_BAC, sd_MADR, sd_BP, sd_HDM, sd_crude,
      sd_full, sd_target, sd_full2, sd_target2,
      included_GBCEEv2, included_GBCEEXv2,
      included_OAL, included_CTMLE, included_GLIDER, included_BAC, included_MADRx, included_MADRy, included_BP, 
      CP_GBCEE_TMLEv2,
      CP_OAL, CP_CTMLE, CP_GLIDER, CP_BAC, CP_MADR, CP_BP, CP_HDM, CP_crude,
      CP_full, CP_target, CP_full2, CP_target2, coefs_OAL_TMLE, sd_OAL_TMLE, CP_OAL_TMLE,
      coefs_SS, sd_SS, included_SS, CP_SS,
      file = file.name);
}
















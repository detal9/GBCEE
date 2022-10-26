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
# BayesPen needs to be downloaded here: https://cran.r-project.org/src/contrib/Archive/BayesPen/
#  and then manually installed
#install.packages("penalized");

### Load required packages;
require(MASS);
require(mvtnorm);
require(BMA);
require(ctmle);
require(boot);
require(bacr);
require(madr);
require(BayesPen);
require(hdm);
require(BCEE);
require(penalized);
source("OAL_V3.R");
source("glider.R");



### Change working directory for saving simulation results
# setwd("...");

############### Scenario 1 - Binary X -- Continuous Y

p = 40;
trueATE = trueB = 1;
var.list = c(paste("U",1:p,sep=""));
nrep = 1;

coefs_GBCEE_TMLEv2 = numeric(nrep);
coefs_OAL = numeric(nrep);
coefs_OAL_TMLE = numeric(nrep);
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

sd_GBCEE_TMLEv2 = numeric(nrep);
sd_OAL = numeric(nrep);
sd_OAL_TMLE = numeric(nrep);
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

included_GBCEEv2 = matrix(, nrow = nrep, ncol = p);
included_GBCEEXv2 = matrix(, nrow = nrep, ncol = p);
included_OAL = matrix(, nrow = nrep, ncol = p);
included_CTMLE = matrix(, nrow = nrep, ncol = p);
included_GLIDER = matrix(, nrow = nrep, ncol = p);
included_BAC = matrix(, nrow = nrep, ncol = p);
included_MADRx = matrix(, nrow = nrep, ncol = p);
included_MADRy = matrix(, nrow = nrep, ncol = p);
included_BP = matrix(, nrow = nrep, ncol = p);

CP_GBCEE_TMLEv2 = numeric(nrep);
CP_OAL = numeric(nrep);
CP_OAL_TMLE = numeric(nrep);
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

stat.madr = function(x, i)
{
 Y = as.matrix(ds$Y[i]);
 X = as.matrix(ds$X[i]);
 U = as.matrix(ds[i,-c(1,2)]);
 madr(Y = Y, X = X, U = U)$madr;
}

for(n in 1000)
{
 set.seed(61947919);
 print(n);
 for(i in 1:nrep)
 {

  ### simulate data
  U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1, nrow = p) + matrix(0, nrow = p, ncol = p));
  U[,1:5] = 1*rowSums(U[,11:15]) + matrix(rnorm(n*5), ncol = 5);
  X = rbinom(n, size = 1, p = plogis(1*rowSums(U[,c(11:30)])));
  Y = rnorm(n, mean = 0.1*rowSums(U[,1:10]) + trueB*X, sd = 1);

  a = Sys.time();
  ### BCEE (omega = 500*n**0.5)
  results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20);
  coefs_GBCEE_TMLEv2[i] = results$beta[1];
  sd_GBCEE_TMLEv2[i] = results$stderr[1];
  included_GBCEEv2[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  included_GBCEEXv2[i, ] = colSums(results$models.X[,1:p]*results$models.X[,p+1]);
  CP_GBCEE_TMLEv2[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];
  b = Sys.time();
  cat("GBCEE EIF variance execution time was ", b-a, "seconds\n");

  a = Sys.time();
  ### BCEE (omega = 500*n**0.5)
  results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "binomial",
		niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "gaussian", X0 = 0, X1 = 1, OR = 20,
		var.comp = "bootstrap");
  coefs_GBCEE_TMLEv2[i] = results$beta[1];
  sd_GBCEE_TMLEv2[i] = results$stderr[1];
  included_GBCEEv2[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+3]);
  included_GBCEEXv2[i, ] = colSums(results$models.X[,1:p]*results$models.X[,p+1]);
  CP_GBCEE_TMLEv2[i] = trueATE > results$beta[1] - 1.96*results$stderr[1] & trueATE < results$beta[1] + 1.96*results$stderr[1];
  b = Sys.time();
  difft = difftime(b, a, units = "secs");
  cat("GBCEE bootstrap variance execution time was ", difft, "seconds\n");

  ### Outcome-adaptive lasso
  a = Sys.time();
  results = OAL.lm(Y = Y, A = X, X = U);
  b = Sys.time();
  coefs_OAL[i] = results$ATE;
  included_OAL[i,] = results$included;
  d = Sys.time();
  tryCatch({
  results = OAL.lm.boot(Y = Y, A = X, X = U, R = 1000);
  CP_OAL[i] = trueATE > quantile(results$t, 0.025) & trueATE < quantile(results$t, 0.975);},
   error = function(e){CP_OAL[i] = NA;});
  e = Sys.time();
  difft1 = difftime(b, a, units = "secs");
  difft2 = difftime(e, d, units = "secs");
  cat("OAL point estimate execution time was ", difft1, "seconds\n",
      "OAL confidence intervals execution time was ", difft2, "seconds\n");


  ### C-TMLE
  a = Sys.time();
  tryCatch({
  results = ctmleDiscrete(Y = Y, A = X, W = data.frame(U), preOrder = FALSE);
  coefs_CTMLE[i] = results$est;
  sd_CTMLE[i] = sqrt(results$var.psi);
  included_CTMLE[i, ] = paste("X", seq(1:p), sep = "") %in% head(results$candidates$terms, results$best_k);
  CP_CTMLE[i] = trueATE >  coefs_CTMLE[i] - 1.96*sd_CTMLE[i] & trueATE <  coefs_CTMLE[i] + 1.96*sd_CTMLE[i];},
   error = function(e){coefs_CTMLE[i] = NA; sd_CTMLE[i] = NA; CP_CTMLE[i]  = NA;}); # included_CTMLE = NA by default
  b = Sys.time();
  difft1 = difftime(b, a, units = "secs");
  cat("C-TMLE execution time was ", difft1, "seconds\n");

  ### GLIDER
  a = Sys.time();
  results = GLiDeR(Xorig = U, Yorig = Y, Aorig = X);
  coefs_GLIDER[i] = results$delta;
  included_GLIDER[i, ] = results$alpha > 0.00001;
  b = Sys.time();
  tryCatch({
  results = GLiDeR_bootstrap(Xorig = U, Yorig = Y, Aorig = X, lambda_star = results$lambda_star, B = 1000);
  CP_GLIDER[i] = trueATE > quantile(results, 0.025) & trueATE < quantile(results, 0.975);},
   error = function(e){CP_GLIDER[i] = NA;});
  d = Sys.time();
  difft1 = difftime(b, a, units = "secs");
  difft2 = difftime(d, a, units = "secs");
  cat("GLiDeR point estimate execution time was ", difft1, "seconds\n",
      "GLiDeR confidence intervals execution time was ", difft2, "seconds\n");


  ### BAC
  ds = data.frame(Y,X,U);
  colnames(ds) = c("Y", "X", paste("U", 1:p, sep = ""));
  varlist = paste("U", 1:p, sep = "");
  a = Sys.time();
  results = bac(data = ds, exposure = "X", outcome = "Y", confounders = varlist, interactors = NULL,
                familyX = "binomial", family = "gaussian", num_its = 5000, burnM = 500, burnB = 500, thin = 5);
  b = Sys.time();
  coefs_BAC[i] = mean(results$ACE);
  sd_BAC[i] = sd(results$ACE);
  included_BAC[i,] = summary(results)$PIP;
  CP_BAC[i] = trueATE < quantile(results$ACE, 0.975) & trueATE > quantile(results$ACE, 0.025);
  difft1 = difftime(b, a, units = "secs");
  cat("BAC execution time was ", difft1, "seconds\n");


  ### MADR
  a = Sys.time();
  capture.output({
  results = madr(Y = Y, X = X, U = U);
  b = Sys.time();
  coefs_MADR[i] = results$madr;
  included_MADRx[i, ] = results$weight.ps;
  included_MADRy[i, ] = results$weight.om;
  d = Sys.time();
  tryCatch({
  results = boot(data.frame(Y,X,U), statistic = stat.madr, R = 1000)$t;
  CP_MADR[i] = trueATE > quantile(results, 0.025) & trueATE < quantile(results, 0.975);},
   error = function(e){CP_MADR[i] = NA;});}, file = 'NUL');
  e = Sys.time();
  difft1 = difftime(b, a, units = "secs");
  difft2 = difftime(e, d, units = "secs");
  cat("MADR point estimate execution time was ", difft1, "seconds\n",
      "MADR confidence intervals execution time was ", difft2, "seconds\n");


  ### BayesPen
  a = Sys.time();
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
  included_BP[i, ] = paste("U", seq(1:p), sep = "") %in% varlist;},
   error = function(e){coefs_BP[i] = NA; sd_BP[i] = NA; CP_BP[i] = NA; included_BP[i, ] = NA;});
  b = Sys.time();
  difft1 = difftime(b, a, units = "secs");
  cat("BP execution time was ", difft1, "seconds\n");

  ### HDM
  a = Sys.time();
  results = rlassoEffect(x = U, y = Y, d = X);
  coefs_HDM[i] = results$coefficients;
  sd_HDM[i] = results$se;
  CP_HDM[i] = trueATE >  coefs_HDM[i] - 1.96*sd_HDM[i] & trueATE <  coefs_HDM[i] + 1.96*sd_HDM[i];
  b = Sys.time();
  difft1 = difftime(b, a, units = "secs");
  cat("HDM execution time was ", difft1, "seconds\n");
 }
}
















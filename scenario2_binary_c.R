expit = plogis;

### Set working directory
# setwd("...");

### Install required packages if necessary
#install.packages("BCEE");

### Load required packages;
require(BCEE);


### Change working directory for saving simulation results
# setwd("...");

### Initializing objects

p = 20;
rho = 0.5;
trueB = 2;
beta = matrix(c(rep(0.6, 4), 0, 0,  rep(0, p - 6)), ncol = 1);
nu = matrix(c(1, 1, 0, 0, 1, 1, rep(0, p - 6)), ncol = 1)*0.6;
var.list = c(paste("U",1:p,sep=""));
nrep = 1000;

coefs_GBCEE = numeric(nrep);
coefs_crude = numeric(nrep);
coefs_full = numeric(nrep);
coefs_target = numeric(nrep);
coefs_full2 = numeric(nrep);
coefs_target2 = numeric(nrep);

sd_GBCEE = numeric(nrep);
sd_crude = numeric(nrep);
sd_full = numeric(nrep);
sd_target = numeric(nrep);
sd_full2 = numeric(nrep);
sd_target2 = numeric(nrep);

included_GBCEE = matrix(, nrow = nrep, ncol = p);
included_GBCEEX = matrix(, nrow = nrep, ncol = p);

CP_GBCEE = numeric(nrep);
CP_crude = numeric(nrep);
CP_full = numeric(nrep);
CP_target = numeric(nrep);
CP_full2 = numeric(nrep);
CP_target2 = numeric(nrep);



n = 1000000;
U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
X = rnorm(n, mean = (U%*%nu));
Xr = rnorm(n, mean = mean(X), sd = sd(X));
YX = rbinom(n, size = 1, p = plogis(trueB*Xr + U%*%beta));
Wmod = glm(YX~Xr, family = binomial(link = "logit"));
True = coef(Wmod)[2];
tronc = 0.99


### Simulation

pb = txtProgressBar(min = 1, max = nrep, initial = 1, style = 3, file = "", char = "*");
n = 1000;
set.seed(61947919);
for(i in 1:nrep)
{
 setTxtProgressBar(pb, i);

 ### simulate data
 U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
 X = rnorm(n, mean = U%*%nu);
 Y = rbinom(n, size = 1, p = plogis(trueB*X + U%*%beta));

 ### BCEE (omega = 500*n**0.5)
 results = GBCEE(X = X, Y = Y, U = U, omega = 500*n**0.5, family.X = "gaussian",
	niter = 2000, priorX = rep(0.5, p), priorY = rep(0.5, p), family.Y = "binomial", OR = 20);
 coefs_GBCEE[i] = results$beta[2];
 sd_GBCEE[i] = results$stderr[2];
 included_GBCEE[i, ] = colSums(results$models.Y[,1:p]*results$models.Y[,p+5]);
 included_GBCEEX[i, ] = colSums(results$models.X[,1:p]*results$models.X[,p+1]);
 CP_GBCEE[i] = True > results$beta[2] - 1.96*results$stderr[2] & True < results$beta[2] + 1.96*results$stderr[2];

 modX = glm(X ~ 1, family = "gaussian");


 ### G-formula - target
 nsamp = 30;
 V = U[,1:4];
 modY = glm(Y ~ X + V, family = binomial(link = "logit"));
 Xs = rnorm(nsamp*n, mean = coef(modX, type = "res"), sd = sd(residuals(modX)));
 Vs = matrix(rep(V, each = nsamp), nrow = n*nsamp, ncol = ncol(V), byrow = FALSE);
 Ys = expit(cbind(1, Xs, Vs)%*%coef(modY));
 coefs_target[i] = coef(glm(Ys ~ Xs, family = binomial(link = "logit")))[2];

 ### G-formula - full
 nsamp = 30;
 V = U;
 modY = glm(Y ~ X + V, family = binomial(link = "logit"));
 Xs = rnorm(nsamp*n, mean = coef(modX, type = "res"), sd = sd(residuals(modX)));
 Vs = matrix(rep(V, each = nsamp), nrow = n*nsamp, ncol = ncol(V), byrow = FALSE);
 Ys = expit(cbind(1, Xs, Vs)%*%coef(modY));
 coefs_full[i] = coef(glm(Ys ~ Xs, family = binomial(link = "logit")))[2];

 logit = qlogis;
 ### Target model - TMLE
 nsamp = 30;
 V = U[,1:4];
 modXU = glm(X ~ V, family = "gaussian");
 modY = glm(Y ~ X + V, family = binomial(link = "logit"));
 Xs = rnorm(nsamp*n, mean = coef(modX, type = "res"), sd = sd(residuals(modX)));
 Vs = matrix(rep(V, each = nsamp), nrow = n*nsamp, ncol = ncol(V), byrow = FALSE);
 Ys = expit(cbind(1, Xs, Vs)%*%coef(modY));
 hX = dt((X - coef(modX))/sd(residuals(modX)), df = 100);
 pX = dt((X - predict(modXU))/sd(residuals(modXU)), df = 100);
 Yp  = predict(modY, type = "res");
 w = hX/pX;
 w = w/mean(w);
 w = pmin(w, quantile(w, 0.99));
 C1A.1 = rep(1, n);
 C1A.2 = X;
 mod.f = glm(Y ~ offset(logit(Yp)) + C1A.1 + C1A.2 - 1, weights = w, family = binomial(link = "logit"));
 Yps = expit(coef(mod.f)[1] + coef(mod.f)[2]*Xs + logit(Ys));
 Yp2 = expit(coef(mod.f)[1] + coef(mod.f)[2]*X + logit(Yp));
 psi.T = coef(glm(Yps ~ Xs, family = binomial(link = "logit")));
 coefs_target2[i] = psi.T[2];
 M = matrix(, nrow = 2, ncol = 2);
 m_1m = expit(cbind(1, Xs)%*%psi.T)*(1 - expit(cbind(1, Xs)%*%psi.T));
 M[1,1] = - mean(m_1m*1*1);
 M[1,2] = - mean(m_1m*1*Xs);
 M[2,1] = - mean(m_1m*Xs*1);
 M[2,2] = - mean(m_1m*Xs*Xs);
 EIC0 = matrix(0, nrow = 2, ncol = n);
 Ys2 = matrix(Ys, nrow = n, ncol = nsamp);
 Yps2 = matrix(Yps, nrow = n, ncol = nsamp);
 Xs2 = matrix(Xs, nrow = n, ncol = nsamp);
 EIC0[1,] = hX/pX*(Y - Yp2)*1 + rowMeans((Ys2 - Yps2)*1); 
 EIC0[2,] = hX/pX*(Y - Yp2)*X + rowMeans((Ys2 - Yps2)*Xs2); 
 EIC = solve(M)%*%EIC0;
 sd_target2[i] = sqrt(diag(var(t(EIC))/n))[2];
 CP_target2[i] = True > coefs_target2[i] - 1.96*sd_target2[i] &
                 True < coefs_target2[i] + 1.96*sd_target2[i];


 ### Full model - TMLE
 nsamp = 30;
 V = U;
 modX = glm(X ~ 1, family = "gaussian");
 modXU = glm(X ~ V, family = "gaussian");
 modY = glm(Y ~ X + V, family = binomial(link = "logit"));
 Xs = rnorm(nsamp*n, mean = coef(modX, type = "res"), sd = sd(residuals(modX)));
 Vs = matrix(rep(V, each = nsamp), nrow = n*nsamp, ncol = ncol(V), byrow = FALSE);
 Ys = expit(cbind(1, Xs, Vs)%*%coef(modY));
 hX = dt((X - coef(modX))/sd(residuals(modX)), df = 100);
 pX = dt((X - predict(modXU))/sd(residuals(modXU)), df = 100);
 Yp  = predict(modY, type = "res");
 w = hX/pX;
 w = w/mean(w);
 w = pmin(w, quantile(w, 0.99));
 C1A.1 = rep(1, n);
 C1A.2 = X;
 mod.f = glm(Y ~ offset(logit(Yp)) + C1A.1 + C1A.2 - 1, weights = w, family = binomial(link = "logit"));
 Yps = expit(coef(mod.f)[1] + coef(mod.f)[2]*Xs + logit(Ys));
 Yp2 = expit(coef(mod.f)[1] + coef(mod.f)[2]*X + logit(Yp));
 psi.T = coef(glm(Yps ~ Xs, family = binomial(link = "logit")));
 coefs_full2[i] = psi.T[2];
 M = matrix(, nrow = 2, ncol = 2);
 m_1m = expit(cbind(1, Xs)%*%psi.T)*(1 - expit(cbind(1, Xs)%*%psi.T));
 M[1,1] = - mean(m_1m*1*1);
 M[1,2] = - mean(m_1m*1*Xs);
 M[2,1] = - mean(m_1m*Xs*1);
 M[2,2] = - mean(m_1m*Xs*Xs);
 EIC0 = matrix(0, nrow = 2, ncol = n);
 Ys2 = matrix(Ys, nrow = n, ncol = nsamp);
 Yps2 = matrix(Yps, nrow = n, ncol = nsamp);
 Xs2 = matrix(Xs, nrow = n, ncol = nsamp);
 EIC0[1,] = hX/pX*(Y - Yp2)*1 + rowMeans((Ys2 - Yps2)*1); 
 EIC0[2,] = hX/pX*(Y - Yp2)*X + rowMeans((Ys2 - Yps2)*Xs2); 
 EIC = solve(M)%*%EIC0;
 sd_full2[i] = sqrt(diag(var(t(EIC))/n))[2];
 CP_full2[i] = True > coefs_full2[i] - 1.96*sd_full2[i] &
               True < coefs_full2[i] + 1.96*sd_full2[i];
 print(i);
 print(Sys.time());
}

 
### Record data
 file.name = paste("results_scenario2binaryC_n1000_tronc","_",gsub("-","",Sys.Date()),".Rdata", sep="");
 save(coefs_GBCEE, coefs_full, coefs_target, coefs_full2, coefs_target2,
     sd_GBCEE, sd_full, sd_target, sd_full2, sd_target2,
     CP_GBCEE, CP_full, CP_target, CP_full2, CP_target2,
     included_GBCEE, included_GBCEEX,
     file = file.name);












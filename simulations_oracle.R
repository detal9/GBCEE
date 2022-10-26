######### Simulations Oracle GBCEE

require(BCEE);

#### Scenario 1 - 1 covariate

set.seed(46131);
nrep = 1000;
ns = c(20, 50, 100, 200, 1000);
prior.oracle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));

post.oracle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));

est.oracle.GBCEE = est.mle.GBCEE = est.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));
sd.oracle.GBCEE = sd.mle.GBCEE = sd.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));
cov.oracle.GBCEE = cov.mle.GBCEE = cov.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));

## True data generating equations
n = 10000000;
U1 = rnorm(n);
X = rbinom(n, 1, plogis(U1));
Y = rnorm(n, X + 0.1*U1, 1);
alphaY1 = lm(Y ~ X);
alphaY2 = lm(Y ~ X + U1);

## True parameters
delta_1.alpha1 = 0.1;
delta_1.alpha2 = 0.1;
sigma_U1 = 1;
sigma_Y = sd(Y);


j = 1;
for(n in ns){
  for(i in 1:nrep){
    U1 = rnorm(n);
    X = rbinom(n, 1, plogis(U1));
    Y = rnorm(n, X + 0.1*U1, 1);

    # Outcome models
    alphaY1 = lm(Y ~ X);
    alphaY2 = lm(Y ~ X + U1);
    BIC.Y1 = BIC(alphaY1);
    BIC.Y2 = BIC(alphaY2);
    PY_alphaY1 = exp(-0.5*BIC.Y1 + n);
    PY_alphaY2 = exp(-0.5*BIC.Y2 + n);

    # Compute posterior of alphaX
    alphaX1 = lm(X ~ 1);
    alphaX2 = lm(X ~ U1);
    BIC.X1 = BIC(alphaX1);
    BIC.X2 = BIC(alphaX2);
    PalphaX1 = exp(-0.5*BIC.X1)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2));
    PalphaX2 = exp(-0.5*BIC.X2)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2));

    # Compute oracle prior
    omega_1.alpha1 = 500*sqrt(n)*(delta_1.alpha1*sigma_U1/sigma_Y)**2;
    omega_1.alpha2 = 500*sqrt(n)*(delta_1.alpha2*sigma_U1/sigma_Y)**2;
    UPalphaY1_alphaX1 = 0.5;
    UPalphaY2_alphaX1 = 0.5;
    UPalphaY1_alphaX2 = 1/(omega_1.alpha1);
    UPalphaY2_alphaX2 = omega_1.alpha2/(1 + omega_1.alpha2);
    PalphaY1_alphaX1 = UPalphaY1_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1);
    PalphaY2_alphaX1 = UPalphaY2_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1);
    PalphaY1_alphaX2 = UPalphaY1_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2);
    PalphaY2_alphaX2 = UPalphaY2_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2);
    prior.oracle.GBCEE.1[i, j] = PalphaY1 = PalphaY1_alphaX1*PalphaX1 + PalphaY1_alphaX2*PalphaX2;
    prior.oracle.GBCEE.2[i, j] = PalphaY2 = PalphaY2_alphaX1*PalphaX1 + PalphaY2_alphaX2*PalphaX2;

    # Compute oracle posterior
    UPostPY_alphaY1 = PY_alphaY1*PalphaY1;
    UPostPY_alphaY2 = PY_alphaY2*PalphaY2;
    post.oracle.GBCEE.1[i, j] = PostPY_alphaY1 = UPostPY_alphaY1/(UPostPY_alphaY1 + UPostPY_alphaY2);
    post.oracle.GBCEE.2[i, j] = PostPY_alphaY2 = UPostPY_alphaY2/(UPostPY_alphaY1 + UPostPY_alphaY2);

    # Compute oracle estimate, sd and CI
    est.oracle.GBCEE[i, j] = coef(alphaY1)[2]*PostPY_alphaY1 + coef(alphaY2)[2]*PostPY_alphaY2;
    sd.oracle.GBCEE[i, j] = sqrt((summary(alphaY1)$coef[2,2]**2 + coef(alphaY1)[2]**2)*PostPY_alphaY1 +
                                 (summary(alphaY2)$coef[2,2]**2 + coef(alphaY2)[2]**2)*PostPY_alphaY2 -
                                 est.oracle.GBCEE[i, j]**2);
    cov.oracle.GBCEE[i, j] = est.oracle.GBCEE[i, j] - 1.96*sd.oracle.GBCEE[i, j] < 1 &
                             est.oracle.GBCEE[i, j] + 1.96*sd.oracle.GBCEE[i, j] > 1;


    # Compute mle prior
    norms = rnorm(1000);
    delta11 = coef(alphaY2)[3] + summary(alphaY2)$coef[3,2]*norms;
    delta12 = coef(alphaY2)[3] + summary(alphaY2)$coef[3,2]*norms;
    omega_1.alpha1 = 500*sqrt(n)*(delta11*sigma_U1/sigma_Y)**2;
    omega_1.alpha2 = 500*sqrt(n)*(delta12*sigma_U1/sigma_Y)**2; 
    UPalphaY1_alphaX1 = 0.5;
    UPalphaY2_alphaX1 = 0.5;
    UPalphaY1_alphaX2 = mean(1/(omega_1.alpha1 + 1));
    UPalphaY2_alphaX2 = mean(omega_1.alpha2/(1 + omega_1.alpha2));
    PalphaY1_alphaX1 = UPalphaY1_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1);
    PalphaY2_alphaX1 = UPalphaY2_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1);
    PalphaY1_alphaX2 = UPalphaY1_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2);
    PalphaY2_alphaX2 = UPalphaY2_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2);
    prior.mle.GBCEE.1[i, j] = PalphaY1 = PalphaY1_alphaX1*PalphaX1 + PalphaY1_alphaX2*PalphaX2;
    prior.mle.GBCEE.2[i, j] = PalphaY2 = PalphaY2_alphaX1*PalphaX1 + PalphaY2_alphaX2*PalphaX2;

    # Compute mle posterior
    UPostPY_alphaY1 = PY_alphaY1*PalphaY1;
    UPostPY_alphaY2 = PY_alphaY2*PalphaY2;
    post.mle.GBCEE.1[i, j] = PostPY_alphaY1 = UPostPY_alphaY1/(UPostPY_alphaY1 + UPostPY_alphaY2);
    post.mle.GBCEE.2[i, j] = PostPY_alphaY2 = UPostPY_alphaY2/(UPostPY_alphaY1 + UPostPY_alphaY2);

    # Compute mle estimate, sd and CI
    est.mle.GBCEE[i, j] = coef(alphaY1)[2]*PostPY_alphaY1 + coef(alphaY2)[2]*PostPY_alphaY2;
    sd.mle.GBCEE[i, j] = sqrt((summary(alphaY1)$coef[2,2]**2 + coef(alphaY1)[2]**2)*PostPY_alphaY1 +
                              (summary(alphaY2)$coef[2,2]**2 + coef(alphaY2)[2]**2)*PostPY_alphaY2 -
                               est.mle.GBCEE[i, j]**2);
    cov.mle.GBCEE[i, j] = est.mle.GBCEE[i, j] - 1.96*sd.mle.GBCEE[i, j] < 1 &
                          est.mle.GBCEE[i, j] + 1.96*sd.mle.GBCEE[i, j] > 1;
  }
  j = j + 1;
  print(data.frame(j, Sys.time()));
}

par(mfrow = c(2,3));
for(j in 1:5)
{
  plot(prior.oracle.GBCEE.2[,j], prior.mle.GBCEE.2[,j], main = paste("P(alphaY2), n =", ns[j]),
       xlab = "oracle", ylab = "mle", xlim = c(0,1), ylim = c(0,1));
}

par(mfrow = c(2,3));
for(j in 1:5)
{
  plot(post.oracle.GBCEE.2[,j], post.mle.GBCEE.2[,j], main = paste("P(alphaY2|O), n =", ns[j]),
       xlab = "oracle", ylab = "mle", xlim = c(0,1), ylim = c(0,1));
}

boxplot(cbind(post.oracle.GBCEE.2[,1], post.mle.GBCEE.2[,1],
              post.oracle.GBCEE.2[,2], post.mle.GBCEE.2[,2],
              post.oracle.GBCEE.2[,3], post.mle.GBCEE.2[,3],
              post.oracle.GBCEE.2[,4], post.mle.GBCEE.2[,4],
              post.oracle.GBCEE.2[,5], post.mle.GBCEE.2[,5]), xaxt = "n");
axis(1, 1:10, c("Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE"));
axis(1, c(1.5, 3.5, 5.5, 7.5, 9.5), line = 2, c("n = 20", "n = 50", "n = 100", "n = 200", "n = 1000"), tick = FALSE);

colMeans(est.oracle.GBCEE);
colMeans(est.mle.GBCEE);
apply(est.oracle.GBCEE, 2, sd);
apply(est.mle.GBCEE, 2, sd);
colMeans(cov.oracle.GBCEE);
colMeans(cov.mle.GBCEE);




#### Scenario 2 - 2 independent covariates

set.seed(46131);
nrep = 1000;
ns = c(20, 50, 100, 200, 1000);
prior.oracle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));

post.oracle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));


est.oracle.GBCEE = est.mle.GBCEE = est.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));
sd.oracle.GBCEE = sd.mle.GBCEE = sd.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));
cov.oracle.GBCEE = cov.mle.GBCEE = cov.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));

## True data generating equations
n = 10000000;
U1 = rnorm(n);
U2 = rnorm(n);
X = rbinom(n, 1, plogis(U1 + U2));
Y = rnorm(n, X + 0.1*U1 + 0.1*U2, 1);
alphaY1 = lm(Y ~ X);
alphaY2 = lm(Y ~ X + U1);
alphaY3 = lm(Y ~ X + U2);
alphaY4 = lm(Y ~ X + U1 + U2);

## True parameters
delta_1.alpha1 = coef(alphaY2)[3];
delta_2.alpha1 = coef(alphaY3)[3];
delta_1.alpha2 = coef(alphaY2)[3];
delta_2.alpha2 = coef(alphaY4)[4];
delta_1.alpha3 = coef(alphaY4)[3];
delta_2.alpha3 = coef(alphaY3)[3];
delta_1.alpha4 = coef(alphaY4)[3];
delta_2.alpha4 = coef(alphaY4)[4];
sigma_U1 = 1;
sigma_U2 = 1;
sigma_Y = sd(Y);


j = 1;
for(n in ns){
  for(i in 1:nrep){
    U1 = rnorm(n);
    U2 = rnorm(n);
    X = rbinom(n, 1, plogis(U1 + U2));
    Y = rnorm(n, X + 0.1*U1 + 0.1*U2, 1);

    # Outcome models
    alphaY1 = lm(Y ~ X);
    alphaY2 = lm(Y ~ X + U1);
    alphaY3 = lm(Y ~ X + U2);
    alphaY4 = lm(Y ~ X + U1 + U2);
    BIC.Y1 = BIC(alphaY1);
    BIC.Y2 = BIC(alphaY2);
    BIC.Y3 = BIC(alphaY3);
    BIC.Y4 = BIC(alphaY4);
    PY_alphaY1 = exp(-0.5*BIC.Y1 + n);
    PY_alphaY2 = exp(-0.5*BIC.Y2 + n);
    PY_alphaY3 = exp(-0.5*BIC.Y3 + n);
    PY_alphaY4 = exp(-0.5*BIC.Y4 + n);

    # Compute posterior of alphaX
    alphaX1 = lm(X ~ 1);
    alphaX2 = lm(X ~ U1);
    alphaX3 = lm(X ~ U2);
    alphaX4 = lm(X ~ U1 + U2);
    BIC.X1 = BIC(alphaX1);
    BIC.X2 = BIC(alphaX2);
    BIC.X3 = BIC(alphaX3);
    BIC.X4 = BIC(alphaX4);
    PalphaX1 = exp(-0.5*BIC.X1)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));
    PalphaX2 = exp(-0.5*BIC.X2)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));
    PalphaX3 = exp(-0.5*BIC.X3)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));
    PalphaX4 = exp(-0.5*BIC.X4)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));

    # Compute oracle prior
    omega_1.alpha1 = 500*sqrt(n)*(delta_1.alpha1*sigma_U1/sigma_Y)**2;
    omega_1.alpha2 = 500*sqrt(n)*(delta_1.alpha2*sigma_U1/sigma_Y)**2;
    omega_1.alpha3 = 500*sqrt(n)*(delta_1.alpha3*sigma_U1/sigma_Y)**2;
    omega_1.alpha4 = 500*sqrt(n)*(delta_1.alpha4*sigma_U1/sigma_Y)**2;
    omega_2.alpha1 = 500*sqrt(n)*(delta_2.alpha1*sigma_U2/sigma_Y)**2;
    omega_2.alpha2 = 500*sqrt(n)*(delta_2.alpha2*sigma_U2/sigma_Y)**2;
    omega_2.alpha3 = 500*sqrt(n)*(delta_2.alpha3*sigma_U2/sigma_Y)**2;
    omega_2.alpha4 = 500*sqrt(n)*(delta_2.alpha4*sigma_U2/sigma_Y)**2;

    UPalphaY1_alphaX1 = 0.25;
    UPalphaY2_alphaX1 = 0.25;
    UPalphaY3_alphaX1 = 0.25;
    UPalphaY4_alphaX1 = 0.25;
    UPalphaY1_alphaX2 = 1/(omega_1.alpha1 + 1)*0.5;
    UPalphaY2_alphaX2 = omega_1.alpha2/(1 + omega_1.alpha2)*0.5;
    UPalphaY3_alphaX2 = 1/(omega_1.alpha3 + 1)*0.5;
    UPalphaY4_alphaX2 = omega_1.alpha4/(1 + omega_1.alpha4)*0.5;
    UPalphaY1_alphaX3 = 0.5*1/(omega_2.alpha1 + 1);
    UPalphaY2_alphaX3 = 0.5*1/(omega_2.alpha2 + 1);
    UPalphaY3_alphaX3 = 0.5*omega_2.alpha3/(1 + omega_2.alpha3);
    UPalphaY4_alphaX3 = 0.5*omega_2.alpha4/(1 + omega_2.alpha4);
    UPalphaY1_alphaX4 = 1/(omega_1.alpha1 + 1)*1/(omega_2.alpha1);
    UPalphaY2_alphaX4 = omega_1.alpha2/(1 + omega_1.alpha2)*1/(omega_2.alpha2 + 1);
    UPalphaY3_alphaX4 = 1/(omega_1.alpha3 + 1)*omega_2.alpha3/(1 + omega_2.alpha3);
    UPalphaY4_alphaX4 = omega_1.alpha4/(1 + omega_1.alpha4)*omega_2.alpha4/(1 + omega_2.alpha4);

    PalphaY1_alphaX1 = UPalphaY1_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY2_alphaX1 = UPalphaY2_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY3_alphaX1 = UPalphaY3_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY4_alphaX1 = UPalphaY4_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY1_alphaX2 = UPalphaY1_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY2_alphaX2 = UPalphaY2_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY3_alphaX2 = UPalphaY3_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY4_alphaX2 = UPalphaY4_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY1_alphaX3 = UPalphaY1_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY2_alphaX3 = UPalphaY2_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY3_alphaX3 = UPalphaY3_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY4_alphaX3 = UPalphaY4_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY1_alphaX4 = UPalphaY1_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY2_alphaX4 = UPalphaY2_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY3_alphaX4 = UPalphaY3_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY4_alphaX4 = UPalphaY4_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);

    prior.oracle.GBCEE.1[i, j] = PalphaY1 = PalphaY1_alphaX1*PalphaX1 + PalphaY1_alphaX2*PalphaX2 + PalphaY1_alphaX3*PalphaX3 + PalphaY1_alphaX4*PalphaX4;
    prior.oracle.GBCEE.2[i, j] = PalphaY2 = PalphaY2_alphaX1*PalphaX1 + PalphaY2_alphaX2*PalphaX2 + PalphaY2_alphaX3*PalphaX3 + PalphaY2_alphaX4*PalphaX4;
    prior.oracle.GBCEE.3[i, j] = PalphaY3 = PalphaY3_alphaX1*PalphaX1 + PalphaY3_alphaX2*PalphaX2 + PalphaY3_alphaX3*PalphaX3 + PalphaY3_alphaX4*PalphaX4;
    prior.oracle.GBCEE.4[i, j] = PalphaY4 = PalphaY4_alphaX1*PalphaX1 + PalphaY4_alphaX2*PalphaX2 + PalphaY4_alphaX3*PalphaX3 + PalphaY4_alphaX4*PalphaX4;


    # Compute oracle posterior
    UPostPY_alphaY1 = PY_alphaY1*PalphaY1;
    UPostPY_alphaY2 = PY_alphaY2*PalphaY2;
    UPostPY_alphaY3 = PY_alphaY3*PalphaY3;
    UPostPY_alphaY4 = PY_alphaY4*PalphaY4;
    post.oracle.GBCEE.1[i, j] = PostPY_alphaY1 = UPostPY_alphaY1/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.oracle.GBCEE.2[i, j] = PostPY_alphaY2 = UPostPY_alphaY2/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.oracle.GBCEE.3[i, j] = PostPY_alphaY3 = UPostPY_alphaY3/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.oracle.GBCEE.4[i, j] = PostPY_alphaY4 = UPostPY_alphaY4/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);


    # Compute oracle estimate, sd and CI
    est.oracle.GBCEE[i, j] = coef(alphaY1)[2]*PostPY_alphaY1 + coef(alphaY2)[2]*PostPY_alphaY2 +
                             coef(alphaY3)[2]*PostPY_alphaY3 + coef(alphaY4)[2]*PostPY_alphaY4;
    sd.oracle.GBCEE[i, j] = sqrt((summary(alphaY1)$coef[2,2]**2 + coef(alphaY1)[2]**2)*PostPY_alphaY1 +
                                 (summary(alphaY2)$coef[2,2]**2 + coef(alphaY2)[2]**2)*PostPY_alphaY2 +
                                 (summary(alphaY3)$coef[2,2]**2 + coef(alphaY3)[2]**2)*PostPY_alphaY3 +
                                 (summary(alphaY4)$coef[2,2]**2 + coef(alphaY4)[2]**2)*PostPY_alphaY4 -
                                 est.oracle.GBCEE[i, j]**2);
    cov.oracle.GBCEE[i, j] = est.oracle.GBCEE[i, j] - 1.96*sd.oracle.GBCEE[i, j] < 1 &
                             est.oracle.GBCEE[i, j] + 1.96*sd.oracle.GBCEE[i, j] > 1;


    # Compute MLE prior
    norms = rnorm(1000);
    delta11 = coef(alphaY2)[3] + summary(alphaY2)$coef[3,2]*norms;
    delta12 = coef(alphaY2)[3] + summary(alphaY2)$coef[3,2]*norms;
    delta13 = coef(alphaY4)[3] + summary(alphaY4)$coef[3,2]*norms;
    delta14 = coef(alphaY4)[3] + summary(alphaY4)$coef[3,2]*norms;
    delta21 = coef(alphaY3)[3] + summary(alphaY3)$coef[3,2]*norms;
    delta22 = coef(alphaY4)[4] + summary(alphaY4)$coef[4,2]*norms;
    delta23 = coef(alphaY3)[3] + summary(alphaY3)$coef[3,2]*norms;
    delta24 = coef(alphaY4)[4] + summary(alphaY4)$coef[4,2]*norms;
    omega_1.alpha1 = 500*sqrt(n)*(delta11*sigma_U1/sigma_Y)**2;
    omega_1.alpha2 = 500*sqrt(n)*(delta12*sigma_U1/sigma_Y)**2;
    omega_1.alpha3 = 500*sqrt(n)*(delta13*sigma_U1/sigma_Y)**2;
    omega_1.alpha4 = 500*sqrt(n)*(delta14*sigma_U1/sigma_Y)**2;
    omega_2.alpha1 = 500*sqrt(n)*(delta21*sigma_U2/sigma_Y)**2;
    omega_2.alpha2 = 500*sqrt(n)*(delta22*sigma_U2/sigma_Y)**2;
    omega_2.alpha3 = 500*sqrt(n)*(delta23*sigma_U2/sigma_Y)**2;
    omega_2.alpha4 = 500*sqrt(n)*(delta24*sigma_U2/sigma_Y)**2;
    UPalphaY1_alphaX1 = 0.25;
    UPalphaY2_alphaX1 = 0.25;
    UPalphaY3_alphaX1 = 0.25;
    UPalphaY4_alphaX1 = 0.25;
    UPalphaY1_alphaX2 = mean(1/(omega_1.alpha1 + 1))*0.5;
    UPalphaY2_alphaX2 = mean(omega_1.alpha2/(1 + omega_1.alpha2))*0.5;
    UPalphaY3_alphaX2 = mean(1/(omega_1.alpha3 + 1))*0.5;
    UPalphaY4_alphaX2 = mean(omega_1.alpha4/(1 + omega_1.alpha4))*0.5;
    UPalphaY1_alphaX3 = 0.5*mean(1/(omega_2.alpha1 + 1));
    UPalphaY2_alphaX3 = 0.5*mean(1/(omega_2.alpha2 + 1));
    UPalphaY3_alphaX3 = 0.5*mean(omega_2.alpha3/(1 + omega_2.alpha3));
    UPalphaY4_alphaX3 = 0.5*mean(omega_2.alpha4/(1 + omega_2.alpha4));
    UPalphaY1_alphaX4 = mean(1/(omega_1.alpha1 + 1))*mean(1/(omega_2.alpha1));
    UPalphaY2_alphaX4 = mean(omega_1.alpha2/(1 + omega_1.alpha2))*mean(1/(omega_2.alpha2 + 1));
    UPalphaY3_alphaX4 = mean(1/(omega_1.alpha3 + 1))*mean(omega_2.alpha3/(1 + omega_2.alpha3));
    UPalphaY4_alphaX4 = mean(omega_1.alpha4/(1 + omega_1.alpha4))*mean(omega_2.alpha4/(1 + omega_2.alpha4));

    PalphaY1_alphaX1 = UPalphaY1_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY2_alphaX1 = UPalphaY2_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY3_alphaX1 = UPalphaY3_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY4_alphaX1 = UPalphaY4_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY1_alphaX2 = UPalphaY1_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY2_alphaX2 = UPalphaY2_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY3_alphaX2 = UPalphaY3_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY4_alphaX2 = UPalphaY4_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY1_alphaX3 = UPalphaY1_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY2_alphaX3 = UPalphaY2_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY3_alphaX3 = UPalphaY3_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY4_alphaX3 = UPalphaY4_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY1_alphaX4 = UPalphaY1_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY2_alphaX4 = UPalphaY2_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY3_alphaX4 = UPalphaY3_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY4_alphaX4 = UPalphaY4_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);

    prior.mle.GBCEE.1[i, j] = PalphaY1 = PalphaY1_alphaX1*PalphaX1 + PalphaY1_alphaX2*PalphaX2 + PalphaY1_alphaX3*PalphaX3 + PalphaY1_alphaX4*PalphaX4;
    prior.mle.GBCEE.2[i, j] = PalphaY2 = PalphaY2_alphaX1*PalphaX1 + PalphaY2_alphaX2*PalphaX2 + PalphaY2_alphaX3*PalphaX3 + PalphaY2_alphaX4*PalphaX4;
    prior.mle.GBCEE.3[i, j] = PalphaY3 = PalphaY3_alphaX1*PalphaX1 + PalphaY3_alphaX2*PalphaX2 + PalphaY3_alphaX3*PalphaX3 + PalphaY3_alphaX4*PalphaX4;
    prior.mle.GBCEE.4[i, j] = PalphaY4 = PalphaY4_alphaX1*PalphaX1 + PalphaY4_alphaX2*PalphaX2 + PalphaY4_alphaX3*PalphaX3 + PalphaY4_alphaX4*PalphaX4;


    # Compute MLE posterior
    UPostPY_alphaY1 = PY_alphaY1*PalphaY1;
    UPostPY_alphaY2 = PY_alphaY2*PalphaY2;
    UPostPY_alphaY3 = PY_alphaY3*PalphaY3;
    UPostPY_alphaY4 = PY_alphaY4*PalphaY4;
    post.mle.GBCEE.1[i, j] = PostPY_alphaY1 = UPostPY_alphaY1/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.mle.GBCEE.2[i, j] = PostPY_alphaY2 = UPostPY_alphaY2/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.mle.GBCEE.3[i, j] = PostPY_alphaY3 = UPostPY_alphaY3/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.mle.GBCEE.4[i, j] = PostPY_alphaY4 = UPostPY_alphaY4/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);


    # Compute MLE estimate, sd and CI
    est.mle.GBCEE[i, j] = coef(alphaY1)[2]*PostPY_alphaY1 + coef(alphaY2)[2]*PostPY_alphaY2 +
                          coef(alphaY3)[2]*PostPY_alphaY3 + coef(alphaY4)[2]*PostPY_alphaY4;
    sd.mle.GBCEE[i, j] = sqrt((summary(alphaY1)$coef[2,2]**2 + coef(alphaY1)[2]**2)*PostPY_alphaY1 +
                              (summary(alphaY2)$coef[2,2]**2 + coef(alphaY2)[2]**2)*PostPY_alphaY2 +
                              (summary(alphaY3)$coef[2,2]**2 + coef(alphaY3)[2]**2)*PostPY_alphaY3 +
                              (summary(alphaY4)$coef[2,2]**2 + coef(alphaY4)[2]**2)*PostPY_alphaY4 -
                              est.mle.GBCEE[i, j]**2);
    cov.mle.GBCEE[i, j] = est.mle.GBCEE[i, j] - 1.96*sd.mle.GBCEE[i, j] < 1 &
                          est.mle.GBCEE[i, j] + 1.96*sd.mle.GBCEE[i, j] > 1;
  }
  j = j + 1;
  print(data.frame(j, Sys.time()));
}

boxplot(cbind(post.oracle.GBCEE.4[,1], post.mle.GBCEE.4[,1],
              post.oracle.GBCEE.4[,2], post.mle.GBCEE.4[,2],
              post.oracle.GBCEE.4[,3], post.mle.GBCEE.4[,3],
              post.oracle.GBCEE.4[,4], post.mle.GBCEE.4[,4],
              post.oracle.GBCEE.4[,5], post.mle.GBCEE.4[,5]), xaxt = "n");
axis(1, 1:10, c("Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE"));
axis(1, c(1.5, 3.5, 5.5, 7.5, 9.5), line = 2, c("n = 20", "n = 50", "n = 100", "n = 200", "n = 1000"), tick = FALSE);

round(colMeans(est.oracle.GBCEE) - 1, 2);
round(colMeans(est.mle.GBCEE) - 1, 2);
round(apply(est.oracle.GBCEE, 2, sd), 2);
round(apply(est.mle.GBCEE, 2, sd), 2);
round(colMeans(cov.oracle.GBCEE), 2);
round(colMeans(cov.mle.GBCEE), 2);





#### Scenario 3 - 2 covariates with conditional instrument

set.seed(46131);
nrep = 1000;
ns = c(20, 50, 100, 200, 1000);
prior.oracle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
prior.oracle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
prior.mle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
prior.approx.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));

post.oracle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
post.oracle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
post.mle.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.1 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.2 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.3 = matrix(, nrow = nrep, ncol = length(ns));
post.approx.GBCEE.4 = matrix(, nrow = nrep, ncol = length(ns));


est.oracle.GBCEE = est.mle.GBCEE = est.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));
sd.oracle.GBCEE = sd.mle.GBCEE = sd.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));
cov.oracle.GBCEE = cov.mle.GBCEE = cov.approx.GBCEE = matrix(, nrow = nrep, ncol = length(ns));

## True data generating equations
n = 10000000;
U1 = rnorm(n);
U2 = U1 + rnorm(n);
X = rbinom(n, 1, plogis(U1));
Y = rnorm(n, X + 0.1*U2, 1);
alphaY1 = lm(Y ~ X);
alphaY2 = lm(Y ~ X + U1);
alphaY3 = lm(Y ~ X + U2);
alphaY4 = lm(Y ~ X + U1 + U2);


## True parameters
delta_1.alpha1 = coef(alphaY2)[3];
delta_2.alpha1 = coef(alphaY3)[3];
delta_1.alpha2 = coef(alphaY2)[3];
delta_2.alpha2 = coef(alphaY4)[4];
delta_1.alpha3 = coef(alphaY4)[3];
delta_2.alpha3 = coef(alphaY3)[3];
delta_1.alpha4 = coef(alphaY4)[3];
delta_2.alpha4 = coef(alphaY4)[4];
sigma_U1 = 1;
sigma_U2 = sd(U2);
sigma_Y = sd(Y);


j = 1;
for(n in ns){
  for(i in 1:nrep){
    U1 = rnorm(n);
    U2 = U1 + rnorm(n);
    X = rbinom(n, 1, plogis(U1));
    Y = rnorm(n, X + 0.1*U2, 1);

    # Outcome models
    alphaY1 = lm(Y ~ X);
    alphaY2 = lm(Y ~ X + U1);
    alphaY3 = lm(Y ~ X + U2);
    alphaY4 = lm(Y ~ X + U1 + U2);
    BIC.Y1 = BIC(alphaY1);
    BIC.Y2 = BIC(alphaY2);
    BIC.Y3 = BIC(alphaY3);
    BIC.Y4 = BIC(alphaY4);
    PY_alphaY1 = exp(-0.5*BIC.Y1 + n);
    PY_alphaY2 = exp(-0.5*BIC.Y2 + n);
    PY_alphaY3 = exp(-0.5*BIC.Y3 + n);
    PY_alphaY4 = exp(-0.5*BIC.Y4 + n);

    # Compute posterior of alphaX
    alphaX1 = lm(X ~ 1);
    alphaX2 = lm(X ~ U1);
    alphaX3 = lm(X ~ U2);
    alphaX4 = lm(X ~ U1 + U2);
    BIC.X1 = BIC(alphaX1);
    BIC.X2 = BIC(alphaX2);
    BIC.X3 = BIC(alphaX3);
    BIC.X4 = BIC(alphaX4);
    PalphaX1 = exp(-0.5*BIC.X1)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));
    PalphaX2 = exp(-0.5*BIC.X2)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));
    PalphaX3 = exp(-0.5*BIC.X3)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));
    PalphaX4 = exp(-0.5*BIC.X4)/(exp(-0.5*BIC.X1) + exp(-0.5*BIC.X2) + exp(-0.5*BIC.X3) + exp(-0.5*BIC.X4));

    # Compute oracle prior
    omega_1.alpha1 = 500*sqrt(n)*(delta_1.alpha1*sigma_U1/sigma_Y)**2;
    omega_1.alpha2 = 500*sqrt(n)*(delta_1.alpha2*sigma_U1/sigma_Y)**2;
    omega_1.alpha3 = 500*sqrt(n)*(delta_1.alpha3*sigma_U1/sigma_Y)**2;
    omega_1.alpha4 = 500*sqrt(n)*(delta_1.alpha4*sigma_U1/sigma_Y)**2;
    omega_2.alpha1 = 500*sqrt(n)*(delta_2.alpha1*sigma_U2/sigma_Y)**2;
    omega_2.alpha2 = 500*sqrt(n)*(delta_2.alpha2*sigma_U2/sigma_Y)**2;
    omega_2.alpha3 = 500*sqrt(n)*(delta_2.alpha3*sigma_U2/sigma_Y)**2;
    omega_2.alpha4 = 500*sqrt(n)*(delta_2.alpha4*sigma_U2/sigma_Y)**2;

    UPalphaY1_alphaX1 = 0.25;
    UPalphaY2_alphaX1 = 0.25;
    UPalphaY3_alphaX1 = 0.25;
    UPalphaY4_alphaX1 = 0.25;
    UPalphaY1_alphaX2 = 1/(omega_1.alpha1 + 1)*0.5;
    UPalphaY2_alphaX2 = omega_1.alpha2/(1 + omega_1.alpha2)*0.5;
    UPalphaY3_alphaX2 = 1/(omega_1.alpha3 + 1)*0.5;
    UPalphaY4_alphaX2 = omega_1.alpha4/(1 + omega_1.alpha4)*0.5;
    UPalphaY1_alphaX3 = 0.5*1/(omega_2.alpha1 + 1);
    UPalphaY2_alphaX3 = 0.5*1/(omega_2.alpha2 + 1);
    UPalphaY3_alphaX3 = 0.5*omega_2.alpha3/(1 + omega_2.alpha3);
    UPalphaY4_alphaX3 = 0.5*omega_2.alpha4/(1 + omega_2.alpha4);
    UPalphaY1_alphaX4 = 1/(omega_1.alpha1 + 1)*1/(omega_2.alpha1);
    UPalphaY2_alphaX4 = omega_1.alpha2/(1 + omega_1.alpha2)*1/(omega_2.alpha2 + 1);
    UPalphaY3_alphaX4 = 1/(omega_1.alpha3 + 1)*omega_2.alpha3/(1 + omega_2.alpha3);
    UPalphaY4_alphaX4 = omega_1.alpha4/(1 + omega_1.alpha4)*omega_2.alpha4/(1 + omega_2.alpha4);

    PalphaY1_alphaX1 = UPalphaY1_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY2_alphaX1 = UPalphaY2_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY3_alphaX1 = UPalphaY3_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY4_alphaX1 = UPalphaY4_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY1_alphaX2 = UPalphaY1_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY2_alphaX2 = UPalphaY2_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY3_alphaX2 = UPalphaY3_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY4_alphaX2 = UPalphaY4_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY1_alphaX3 = UPalphaY1_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY2_alphaX3 = UPalphaY2_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY3_alphaX3 = UPalphaY3_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY4_alphaX3 = UPalphaY4_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY1_alphaX4 = UPalphaY1_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY2_alphaX4 = UPalphaY2_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY3_alphaX4 = UPalphaY3_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY4_alphaX4 = UPalphaY4_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);

    prior.oracle.GBCEE.1[i, j] = PalphaY1 = PalphaY1_alphaX1*PalphaX1 + PalphaY1_alphaX2*PalphaX2 + PalphaY1_alphaX3*PalphaX3 + PalphaY1_alphaX4*PalphaX4;
    prior.oracle.GBCEE.2[i, j] = PalphaY2 = PalphaY2_alphaX1*PalphaX1 + PalphaY2_alphaX2*PalphaX2 + PalphaY2_alphaX3*PalphaX3 + PalphaY2_alphaX4*PalphaX4;
    prior.oracle.GBCEE.3[i, j] = PalphaY3 = PalphaY3_alphaX1*PalphaX1 + PalphaY3_alphaX2*PalphaX2 + PalphaY3_alphaX3*PalphaX3 + PalphaY3_alphaX4*PalphaX4;
    prior.oracle.GBCEE.4[i, j] = PalphaY4 = PalphaY4_alphaX1*PalphaX1 + PalphaY4_alphaX2*PalphaX2 + PalphaY4_alphaX3*PalphaX3 + PalphaY4_alphaX4*PalphaX4;


    # Compute oracle posterior
    UPostPY_alphaY1 = PY_alphaY1*PalphaY1;
    UPostPY_alphaY2 = PY_alphaY2*PalphaY2;
    UPostPY_alphaY3 = PY_alphaY3*PalphaY3;
    UPostPY_alphaY4 = PY_alphaY4*PalphaY4;
    post.oracle.GBCEE.1[i, j] = PostPY_alphaY1 = UPostPY_alphaY1/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.oracle.GBCEE.2[i, j] = PostPY_alphaY2 = UPostPY_alphaY2/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.oracle.GBCEE.3[i, j] = PostPY_alphaY3 = UPostPY_alphaY3/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.oracle.GBCEE.4[i, j] = PostPY_alphaY4 = UPostPY_alphaY4/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);


    # Compute oracle estimate, sd and CI
    est.oracle.GBCEE[i, j] = coef(alphaY1)[2]*PostPY_alphaY1 + coef(alphaY2)[2]*PostPY_alphaY2 +
                             coef(alphaY3)[2]*PostPY_alphaY3 + coef(alphaY4)[2]*PostPY_alphaY4;
    sd.oracle.GBCEE[i, j] = sqrt((summary(alphaY1)$coef[2,2]**2 + coef(alphaY1)[2]**2)*PostPY_alphaY1 +
                                 (summary(alphaY2)$coef[2,2]**2 + coef(alphaY2)[2]**2)*PostPY_alphaY2 +
                                 (summary(alphaY3)$coef[2,2]**2 + coef(alphaY3)[2]**2)*PostPY_alphaY3 +
                                 (summary(alphaY4)$coef[2,2]**2 + coef(alphaY4)[2]**2)*PostPY_alphaY4 -
                                 est.oracle.GBCEE[i, j]**2);
    cov.oracle.GBCEE[i, j] = est.oracle.GBCEE[i, j] - 1.96*sd.oracle.GBCEE[i, j] < 1 &
                             est.oracle.GBCEE[i, j] + 1.96*sd.oracle.GBCEE[i, j] > 1;


    # Compute MLE prior
    norms = rnorm(1000);
    delta11 = coef(alphaY2)[3] + summary(alphaY2)$coef[3,2]*norms;
    delta12 = coef(alphaY2)[3] + summary(alphaY2)$coef[3,2]*norms;
    delta13 = coef(alphaY4)[3] + summary(alphaY4)$coef[3,2]*norms;
    delta14 = coef(alphaY4)[3] + summary(alphaY4)$coef[3,2]*norms;
    delta21 = coef(alphaY3)[3] + summary(alphaY3)$coef[3,2]*norms;
    delta22 = coef(alphaY4)[4] + summary(alphaY4)$coef[4,2]*norms;
    delta23 = coef(alphaY3)[3] + summary(alphaY3)$coef[3,2]*norms;
    delta24 = coef(alphaY4)[4] + summary(alphaY4)$coef[4,2]*norms;
    omega_1.alpha1 = 500*sqrt(n)*(delta11*sigma_U1/sigma_Y)**2;
    omega_1.alpha2 = 500*sqrt(n)*(delta12*sigma_U1/sigma_Y)**2;
    omega_1.alpha3 = 500*sqrt(n)*(delta13*sigma_U1/sigma_Y)**2;
    omega_1.alpha4 = 500*sqrt(n)*(delta14*sigma_U1/sigma_Y)**2;
    omega_2.alpha1 = 500*sqrt(n)*(delta21*sigma_U2/sigma_Y)**2;
    omega_2.alpha2 = 500*sqrt(n)*(delta22*sigma_U2/sigma_Y)**2;
    omega_2.alpha3 = 500*sqrt(n)*(delta23*sigma_U2/sigma_Y)**2;
    omega_2.alpha4 = 500*sqrt(n)*(delta24*sigma_U2/sigma_Y)**2;

    UPalphaY1_alphaX1 = 0.25;
    UPalphaY2_alphaX1 = 0.25;
    UPalphaY3_alphaX1 = 0.25;
    UPalphaY4_alphaX1 = 0.25;
    UPalphaY1_alphaX2 = mean(1/(omega_1.alpha1 + 1))*0.5;
    UPalphaY2_alphaX2 = mean(omega_1.alpha2/(1 + omega_1.alpha2))*0.5;
    UPalphaY3_alphaX2 = mean(1/(omega_1.alpha3 + 1))*0.5;
    UPalphaY4_alphaX2 = mean(omega_1.alpha4/(1 + omega_1.alpha4))*0.5;
    UPalphaY1_alphaX3 = 0.5*mean(1/(omega_2.alpha1 + 1));
    UPalphaY2_alphaX3 = 0.5*mean(1/(omega_2.alpha2 + 1));
    UPalphaY3_alphaX3 = 0.5*mean(omega_2.alpha3/(1 + omega_2.alpha3));
    UPalphaY4_alphaX3 = 0.5*mean(omega_2.alpha4/(1 + omega_2.alpha4));
    UPalphaY1_alphaX4 = mean(1/(omega_1.alpha1 + 1))*mean(1/(omega_2.alpha1));
    UPalphaY2_alphaX4 = mean(omega_1.alpha2/(1 + omega_1.alpha2))*mean(1/(omega_2.alpha2 + 1));
    UPalphaY3_alphaX4 = mean(1/(omega_1.alpha3 + 1))*mean(omega_2.alpha3/(1 + omega_2.alpha3));
    UPalphaY4_alphaX4 = mean(omega_1.alpha4/(1 + omega_1.alpha4))*mean(omega_2.alpha4/(1 + omega_2.alpha4));

    PalphaY1_alphaX1 = UPalphaY1_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY2_alphaX1 = UPalphaY2_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY3_alphaX1 = UPalphaY3_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY4_alphaX1 = UPalphaY4_alphaX1/(UPalphaY1_alphaX1 + UPalphaY2_alphaX1 + UPalphaY3_alphaX1 + UPalphaY4_alphaX1);
    PalphaY1_alphaX2 = UPalphaY1_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY2_alphaX2 = UPalphaY2_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY3_alphaX2 = UPalphaY3_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY4_alphaX2 = UPalphaY4_alphaX2/(UPalphaY1_alphaX2 + UPalphaY2_alphaX2 + UPalphaY3_alphaX2 + UPalphaY4_alphaX2);
    PalphaY1_alphaX3 = UPalphaY1_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY2_alphaX3 = UPalphaY2_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY3_alphaX3 = UPalphaY3_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY4_alphaX3 = UPalphaY4_alphaX2/(UPalphaY1_alphaX3 + UPalphaY2_alphaX3 + UPalphaY3_alphaX3 + UPalphaY4_alphaX3);
    PalphaY1_alphaX4 = UPalphaY1_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY2_alphaX4 = UPalphaY2_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY3_alphaX4 = UPalphaY3_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);
    PalphaY4_alphaX4 = UPalphaY4_alphaX2/(UPalphaY1_alphaX4 + UPalphaY2_alphaX4 + UPalphaY3_alphaX4 + UPalphaY4_alphaX4);

    prior.mle.GBCEE.1[i, j] = PalphaY1 = PalphaY1_alphaX1*PalphaX1 + PalphaY1_alphaX2*PalphaX2 + PalphaY1_alphaX3*PalphaX3 + PalphaY1_alphaX4*PalphaX4;
    prior.mle.GBCEE.2[i, j] = PalphaY2 = PalphaY2_alphaX1*PalphaX1 + PalphaY2_alphaX2*PalphaX2 + PalphaY2_alphaX3*PalphaX3 + PalphaY2_alphaX4*PalphaX4;
    prior.mle.GBCEE.3[i, j] = PalphaY3 = PalphaY3_alphaX1*PalphaX1 + PalphaY3_alphaX2*PalphaX2 + PalphaY3_alphaX3*PalphaX3 + PalphaY3_alphaX4*PalphaX4;
    prior.mle.GBCEE.4[i, j] = PalphaY4 = PalphaY4_alphaX1*PalphaX1 + PalphaY4_alphaX2*PalphaX2 + PalphaY4_alphaX3*PalphaX3 + PalphaY4_alphaX4*PalphaX4;


    # Compute MLE posterior
    UPostPY_alphaY1 = PY_alphaY1*PalphaY1;
    UPostPY_alphaY2 = PY_alphaY2*PalphaY2;
    UPostPY_alphaY3 = PY_alphaY3*PalphaY3;
    UPostPY_alphaY4 = PY_alphaY4*PalphaY4;
    post.mle.GBCEE.1[i, j] = PostPY_alphaY1 = UPostPY_alphaY1/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.mle.GBCEE.2[i, j] = PostPY_alphaY2 = UPostPY_alphaY2/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.mle.GBCEE.3[i, j] = PostPY_alphaY3 = UPostPY_alphaY3/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);
    post.mle.GBCEE.4[i, j] = PostPY_alphaY4 = UPostPY_alphaY4/(UPostPY_alphaY1 + UPostPY_alphaY2 + UPostPY_alphaY3 + UPostPY_alphaY4);


    # Compute MLE estimate, sd and CI
    est.mle.GBCEE[i, j] = coef(alphaY1)[2]*PostPY_alphaY1 + coef(alphaY2)[2]*PostPY_alphaY2 +
                          coef(alphaY3)[2]*PostPY_alphaY3 + coef(alphaY4)[2]*PostPY_alphaY4;
    sd.mle.GBCEE[i, j] = sqrt((summary(alphaY1)$coef[2,2]**2 + coef(alphaY1)[2]**2)*PostPY_alphaY1 +
                              (summary(alphaY2)$coef[2,2]**2 + coef(alphaY2)[2]**2)*PostPY_alphaY2 +
                              (summary(alphaY3)$coef[2,2]**2 + coef(alphaY3)[2]**2)*PostPY_alphaY3 +
                              (summary(alphaY4)$coef[2,2]**2 + coef(alphaY4)[2]**2)*PostPY_alphaY4 -
                              est.mle.GBCEE[i, j]**2);
    cov.mle.GBCEE[i, j] = est.mle.GBCEE[i, j] - 1.96*sd.mle.GBCEE[i, j] < 1 &
                          est.mle.GBCEE[i, j] + 1.96*sd.mle.GBCEE[i, j] > 1;
  }
  j = j + 1;
  print(data.frame(j, Sys.time()));
}

par(mfrow = c(2,2));
boxplot(cbind(post.oracle.GBCEE.1[,1], post.mle.GBCEE.1[,1],
              post.oracle.GBCEE.1[,2], post.mle.GBCEE.1[,2],
              post.oracle.GBCEE.1[,3], post.mle.GBCEE.1[,3],
              post.oracle.GBCEE.1[,4], post.mle.GBCEE.1[,4],
              post.oracle.GBCEE.1[,5], post.mle.GBCEE.1[,5]), xaxt = "n", ylim = c(0,1));
axis(1, 1:10, c("Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE"));
axis(1, c(1.5, 3.5, 5.5, 7.5, 9.5), line = 2, c("n = 20", "n = 50", "n = 100", "n = 200", "n = 1000"), tick = FALSE);

boxplot(cbind(post.oracle.GBCEE.2[,1], post.mle.GBCEE.2[,1],
              post.oracle.GBCEE.2[,2], post.mle.GBCEE.2[,2],
              post.oracle.GBCEE.2[,3], post.mle.GBCEE.2[,3],
              post.oracle.GBCEE.2[,4], post.mle.GBCEE.2[,4],
              post.oracle.GBCEE.2[,5], post.mle.GBCEE.2[,5]), xaxt = "n", ylim = c(0,1));
axis(1, 1:10, c("Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE"));
axis(1, c(1.5, 3.5, 5.5, 7.5, 9.5), line = 2, c("n = 20", "n = 50", "n = 100", "n = 200", "n = 1000"), tick = FALSE);

boxplot(cbind(post.oracle.GBCEE.3[,1], post.mle.GBCEE.3[,1],
              post.oracle.GBCEE.3[,2], post.mle.GBCEE.3[,2],
              post.oracle.GBCEE.3[,3], post.mle.GBCEE.3[,3],
              post.oracle.GBCEE.3[,4], post.mle.GBCEE.3[,4],
              post.oracle.GBCEE.3[,5], post.mle.GBCEE.3[,5]), xaxt = "n", ylim = c(0,1));
axis(1, 1:10, c("Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE"));
axis(1, c(1.5, 3.5, 5.5, 7.5, 9.5), line = 2, c("n = 20", "n = 50", "n = 100", "n = 200", "n = 1000"), tick = FALSE);

boxplot(cbind(post.oracle.GBCEE.4[,1], post.mle.GBCEE.4[,1],
              post.oracle.GBCEE.4[,2], post.mle.GBCEE.4[,2],
              post.oracle.GBCEE.4[,3], post.mle.GBCEE.4[,3],
              post.oracle.GBCEE.4[,4], post.mle.GBCEE.4[,4],
              post.oracle.GBCEE.4[,5], post.mle.GBCEE.4[,5]), xaxt = "n", ylim = c(0,1));
axis(1, 1:10, c("Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE", "Oracle", "MLE"));
axis(1, c(1.5, 3.5, 5.5, 7.5, 9.5), line = 2, c("n = 20", "n = 50", "n = 100", "n = 200", "n = 1000"), tick = FALSE);

round(colMeans(est.oracle.GBCEE) - 1, 2);
round(colMeans(est.mle.GBCEE) - 1, 2);
round(apply(est.oracle.GBCEE, 2, sd), 2);
round(apply(est.mle.GBCEE, 2, sd), 2);
round(colMeans(cov.oracle.GBCEE), 2);
round(colMeans(cov.mle.GBCEE), 2);








require(mvtnorm);
n = 1000000;
layout(mat = matrix(c(1, 3, 5, 7, 9, 11, 13, 
                      2, 4, 6, 8, 10, 12, 14), nrow = 2, ncol = 7, byrow = T));
par(mar = c(2, 4.5, 3, 2))
### Scenario 1
p = 40;
trueATE = trueB = 1;
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1, nrow = p) + matrix(0, nrow = p, ncol = p));
U[,1:5] = 1*rowSums(U[,11:15]) + matrix(rnorm(n*5), ncol = 5);
p1 = plogis(1*rowSums(U[,c(11:30)]));
X = rbinom(n, size = 1, p = p1);
boxplot(p1~X, main = "Scenario 1", ylab = "True exposure model", xlab = "", cex.lab = 2);
boxplot(p1t~X, ylab = "Target exposure model", xlab = "", cex.lab = 2);



### Scenario 2
p = 20;
rho = 0.5;
trueATE = trueB = 2;
beta = matrix(c(rep(0.6, 4), 0, 0,  rep(0, p - 6)), ncol = 1);
nu = matrix(c(1, 1, 0, 0, 1, 1, rep(0, p - 6)), ncol = 1);
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
p1 = plogis(U%*%nu);
X = rbinom(n, size = 1, p = plogis(U%*%nu));
p1t = glm(X~U[,1:4], family = binomial)$fitted;
boxplot(p1~X, main = "Scenario 2", ylab = "", xlab = "");
boxplot(p1t~X, main = "", ylab = "", xlab = "");


### Scenario 3
p = 100;
trueATE = trueB = 1;
beta = matrix(c(2, 0.2, 5, 5, rep(0, p - 4)), ncol = 1);
nu = matrix(c(0.5, -1, 0, 0, 0.3, -0.3, 0.3, -0.3, rep(0, p - 8)), ncol = 1);
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(1, p), sigma = diag(4, nrow = p) + matrix(0, nrow = p, ncol = p));
p1 = plogis(U%*%nu);
X = rbinom(n, size = 1, p = p1);
Y = rnorm(n, mean = trueB*X + U%*%beta, sd = 2);
p1t = glm(X~U[,1:4], family = binomial)$fitted;
boxplot(p1~X, main = "Scenario 3", ylab = "", xlab = "");
boxplot(p1t~X, ylab = "", xlab = "");



### Scenario 4
p = 5;
trueATE = trueB = 1;
rho = 0.6;
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(1, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
U2 = t(apply(U, 1, function(x){c(x%*%t(x))}));
p1 = plogis(0.5*U[,1] + 0.5*U[,2] + 0.1*U[,3]);
X = rbinom(n, size = 1, p = p1);
Y = rnorm(n, mean = trueB*X + U[,3] + U[,4] + U[,5] + rowSums(U2), sd = 1);
p1t = glm(X~U[,1:5], family = binomial)$fitted;
boxplot(p1~X, main = "Scenario 4", ylab = "", xlab = "");
boxplot(p1t~X, ylab = "", xlab = "");


### Scenario 5
p = 5;
trueATE = trueB = 1;
rho = 0.6;
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(1, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
U2 = t(apply(U, 1, function(x){c(x%*%t(x))}));
p1 = plogis(-5 + U[,3] + U[,4] + U[,5] + rowSums(U2));
X = rbinom(n, size = 1, p = plogis(-5 + U[,3] + U[,4] + U[,5] + rowSums(U2)));
Y = rnorm(n, mean = trueB*X + 0.5*U[,1] + 0.5*U[,2] + 0.1*U[,3], sd = 1);
p1t = glm(X ~ U[,3] + U[,4] + U[,5] + 
              U[,3]*(U[,3] + U[,4] + U[,5]) +
              U[,4]*(U[,4] + U[,5]) +
              U[,5]*(U[,5]), family = binomial)$fitted;
boxplot(p1~X, main = "Scenario 5", ylab = "", xlab = "");
boxplot(p1t~X, ylab = "", xlab = "");


### Scenario 6
p = 30;
rho = 0.2;
trueATE = trueB = 2;
beta = matrix(c(rep(0.05, 25), rep(0, p - 25)), ncol = 1);
nu = matrix(c(rep(0, p - 25), rep(0.05, 25)), ncol = 1);
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(0, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
p1 = plogis(U%*%nu);
X = rbinom(n, size = 1, p = p1);
Y = rnorm(n, mean = trueB*X + U%*%beta, sd = 1);
p1t = glm(X ~ U[,1:25], family = binomial)$fitted;
boxplot(p1~X, main = "Scenario 6", ylab = "", xlab = "");
boxplot(p1t~X, ylab = "", xlab = "");


### Scenario 7
p = 5;
trueB = 1;
trueATE = 2;
rho = 0.6;
var.list = c(paste("U",1:p,sep=""));
U = rmvnorm(n = n, mean = rep(1, p), sigma = diag(1-rho, nrow = p) + matrix(rho, nrow = p, ncol = p));
p1 = plogis(0.5*U[,1] + 0.5*U[,2] + 0.1*U[,3]);
X = rbinom(n, size = 1, p = plogis(0.5*U[,1] + 0.5*U[,2] + 0.1*U[,3]));
Y = rnorm(n, mean = trueB*X + U[,3] + U[,4] + U[,5] + X*U[,3], sd = 1);
p1t = glm(X ~ U[,3:5], family = binomial)$fitted;
boxplot(p1~X, main = "Scenario 7", ylab = "", xlab = "");
boxplot(p1t~X, ylab = "", xlab = "");

# set vector of possible lambda's to try for oa-lasso
lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = as.character(lambda_vec)

# lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
gamma_convergence_factor = 2

# get the gamma value for each value in the lambda vector that corresponds to convergence factor
gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
names(gamma_vals) = names(lambda_vec)


### define some functions for generating data, ATE estimates, and the wAMD,
logit <- function(p) {log(p)-log(1-p)}
ATE_est = function(fY,fw,fA){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}
RR_est = function(fY, fw, fA){
  t_RR = fY*fw
  tt_RR = ( ( sum(t_RR[fA==1]) / sum(fw[fA==1]) ) / ( sum(t_RR[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_RR);
} 
create_weights = function(fp,fA){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}
create_mweights = function(fp,fA){
  fw <- pmin(fp,1-fp)/(fp)
  fw[fA==0]  <- pmin(fp[fA==0],1-fp[fA==0])/(1-fp[fA==0])
  return(fw)
}

wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){ 
    this.var = paste("w",varlist[jj],sep="") 
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
    diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec = diff_vec * abs(beta) 
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret) 
}

OAL.lm <- function(Y, A, X, matching.w = FALSE){
  
  data = data.frame(Y, A, X);
  var.list = paste("X", 1:ncol(X), sep = "");
  names(data) = c("Y", "A", var.list);
  n = nrow(data);
  
  # Normalize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(data),byrow=TRUE)
  data[,var.list] = data[,var.list] - Temp.mean
  temp.sd = apply(data[,var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(data),byrow=TRUE)
  data[,var.list] = data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  glm.Y = lm(y.form,data=data);
  betaXY = coef(glm.Y)[var.list] 
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value 
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  #lil <- "-5"
  #lil <- "-2"
  #lil <- "-1"
  #lil <- "-0.75" 
  #lil <- "-0.5"
  #lil <- "-0.25"
  #lil <- "0.25"
  #lil <- "0.49" 
  
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = n^(il)*abs(betaXY)^(-ig);
    ### run outcome-adaptive lasso model with appropriate penalty
    logit_oal = penalized(response = A, penalized = X, data=data, lambda1=oal_pen, model="logistic", trace = FALSE )
    # generate propensity score
    data[,paste("f.pA",lil,sep="")] = fitted(logit_oal)
    # save propensity score coefficients
    coeff_XA[,lil] = (abs(coefficients(logit_oal, "all")[-1]) > 0.0001);
    # create inverse probability of treatment weights
    if(matching.w == FALSE)
    {
     data[,paste("w",lil,sep="")] = create_weights(fp=data[,paste("f.pA",lil,sep="")],fA=data$A)
    } else # matching.w == TRUE
    {
     data[,paste("w",lil,sep="")] = create_mweights(fp=data[,paste("f.pA",lil,sep="")],fA=data$A)
    }
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=data,varlist=var.list,trt.var="A",
                                  wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
    # save ATE estimate for this lambda value
    ATE[lil] = ATE_est(fY=data$Y,fw=data[,paste("w",lil,sep="")],fA=data$A)
    
    
  } # close loop through lambda values
  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)
  W = as.matrix(X[, coeff_XA[,tt], drop = FALSE]);
  TMLE = NA;
  tryCatch({
    if(sum(coeff_XA[,tt]) == 0){
      Q1.0 = mean(Y[A==1]);
      Q0.0 = mean(Y[A==0]);
      Q.0 = mean(Y);
      g1W = mean(A);
      bounds = quantile(g1W, c(0.01, 0.99)); 
      g1W = pmin(g1W, bounds[2]);
      g1W = pmax(g1W, bounds[1]);
      H1 = A/g1W; H0 = (1 - A)/(1 - g1W);
      w = H1 + H0;
      epsilon = coef(lm(Y~1+A+offset(Q.0), weights = w));
      Q1.1 = Q1.0 + sum(epsilon);
      Q0.1 = Q0.0 + epsilon[1];
      TMLE = mean(Q1.1 - Q0.1);
    } else{
      datTMLE = data.frame(A, W);
      Qmod = lm(Y~A+W, data = datTMLE);
      Q1.0 = predict(Qmod, newdata = data.frame(A = 1, W));
      Q0.0 = predict(Qmod, newdata = data.frame(A = 0, W));
      Q.0 = predict(Qmod);
      gmod = glm(A~W, family = "binomial", data = datTMLE);
      g1W = gmod$fitted;
      bounds = quantile(g1W, c(0.01, 0.99)); 
      g1W = pmin(g1W, bounds[2]);
      g1W = pmax(g1W, bounds[1]);
      H1 = A/g1W; H0 = (1 - A)/(1 - g1W);
      w = H1 + H0;
      epsilon = coef(lm(Y~1+A+offset(Q.0), weights = w));
      Q1.1 = Q1.0 + sum(epsilon);
      Q0.1 = Q0.0 + epsilon[1];
      TMLE = mean(Q1.1 - Q0.1);
    }
  }, error = function(e){});
  return(list(ATE = ATE[tt], TMLE = TMLE, included = coeff_XA[,tt]));
}

OAL.logistic <- function(Y, A, X, matching.w = FALSE){
  
  data = data.frame(Y, A, X);
  var.list = paste("X", 1:ncol(X), sep = "");
  colnames(data) = c("Y", "A", var.list);
  
  # Normalize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(data),byrow=TRUE)
  data[,var.list] = data[,var.list] - Temp.mean
  temp.sd = apply(data[,var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(data),byrow=TRUE)
  data[,var.list] = data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  glm.Y = glm(y.form,data=data, family = "binomial");
  betaXY = coef(glm.Y)[var.list] 
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = TMLE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value 
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  #lil <- "-5"
  #lil <- "-2"
  #lil <- "-1"
  #lil <- "-0.75" 
  #lil <- "-0.5"
  #lil <- "-0.25"
  #lil <- "0.25"
  #lil <- "0.49" 
  
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
    oal_pen = n^(il)*abs(betaXY)^(-ig);
    ### run outcome-adaptive lasso model with appropriate penalty
    logit_oal = penalized(response = A, penalized = X, data=data, lambda1=oal_pen, model="logistic", trace = FALSE)
    # generate propensity score
    data[,paste("f.pA",lil,sep="")] = fitted(logit_oal)
    # save propensity score coefficients
    coeff_XA[,lil] = (abs(coefficients(logit_oal, "all")[-1]) > 0.0001);
    # create inverse probability of treatment weights
    if(matching.w == FALSE)
    {
     data[,paste("w",lil,sep="")] = create_weights(fp=data[,paste("f.pA",lil,sep="")],fA=data$A)
    } else # matching.w == TRUE
    {
     data[,paste("w",lil,sep="")] = create_mweights(fp=data[,paste("f.pA",lil,sep="")],fA=data$A)
    }
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=data,varlist=var.list,trt.var="A",
                                  wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
    # save ATE estimate for this lambda value
    ATE[lil] = ATE_est(fY=data$Y,fw=data[,paste("w",lil,sep="")],fA=data$A)
    
    
  } # close loop through lambda values
  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)
  W = X[, coeff_XA[,tt], drop = FALSE];
  TMLE = NA;
  tryCatch({
    if(sum(coeff_XA[,tt]) == 0){
      Q1.0 = mean(Y[A==1]);
      Q0.0 = mean(Y[A==0]);
      Q.0 = mean(Y);
      g1W = mean(A);
      bounds = quantile(g1W, c(0.01, 0.99)); 
      g1W = pmin(g1W, bounds[2]);
      g1W = pmax(g1W, bounds[1]);
      H1 = A/g1W; H0 = (1 - A)/(1 - g1W);
      w = H1 + H0;
      epsilon = coef(glm(Y~1+A+offset(qlogis(Q.0)), weights = w, family = "binomial"));
      Q1.1 = plogis(qlogis(Q1.0) + sum(epsilon));
      Q0.1 = plogis(qlogis(Q0.0) + epsilon[1]);
      TMLE = mean(Q1.1 - Q0.1);
    } else{
      datTMLE = data.frame(A, W);
      Qmod = glm(Y~A+W, data = datTMLE, family = "binomial");
      Q1.0 = predict(Qmod, newdata = data.frame(A = 1, W), type = "res");
      Q0.0 = predict(Qmod, newdata = data.frame(A = 0, W), type = "res");
      Q.0 = predict(Qmod, type = "res");
      gmod = glm(A~W, family = "binomial", data = datTMLE);
      g1W = gmod$fitted;
      bounds = quantile(g1W, c(0.01, 0.99)); 
      g1W = pmin(g1W, bounds[2]);
      g1W = pmax(g1W, bounds[1]);
      H1 = A/g1W; H0 = (1 - A)/(1 - g1W);
      w = H1 + H0;
      epsilon = coef(glm(Y~1+A+offset(qlogis(Q.0)), weights = w, family = "binomial"));
      Q1.1 = plogis(qlogis(Q1.0) + sum(epsilon));
      Q0.1 = plogis(qlogis(Q0.0) + epsilon[1]);
      TMLE = mean(Q1.1 - Q0.1);
    }
  }, error = function(e){});

  return(list(ATE = ATE[tt], TMLE = TMLE, included = coeff_XA[,tt]));
}

OAL.lm.boot = function(Y, A, X, R, matching.w = FALSE)
{
 ds = data.frame(Y, A, X);
 var.list = paste("X", 1:ncol(X), sep = "");
 colnames(ds) = c("Y", "A", var.list);
 bootf = function(ds, i){
   res = OAL.lm(Y = ds$Y[i], A = ds$A[i], X = ds[i,-c(1,2)], matching.w);
   return(c(res$ATE, res$TMLE));
 }
 return(boot(ds, bootf, R));
}

OAL.logistic.boot = function(Y, A, X, R, matching.w = FALSE)
{
 ds = data.frame(Y, A, X);
 var.list = paste("X", 1:ncol(X), sep = "");
 colnames(ds) = c("Y", "A", var.list);
 bootf = function(ds, i){
   res = OAL.logistic(Y = ds$Y[i], A = ds$A[i], X = ds[i,-c(1,2)], matching.w);
   return(c(res$ATE, res$TMLE));
 }
 return(boot(ds, bootf, R));
}

smoothed.boot = function(data, statistic, R)
{
 n = nrow(data);
 index = matrix(sample(n, size = n*R, replace = TRUE), nrow = R, ncol = n);
 stat = apply(index, 1, FUN = function(x,y){statistic(y,x)}, data);
 m.stat = mean(stat);
 yij = apply(index, 1, FUN = function(x,y){table(factor(x, levels = 1:y))}, n); 
 yj = rowMeans(yij);
 yj.mat = matrix(yj, nrow = n, ncol = R, byrow = F);
 stat.mat = matrix(stat, nrow = n, ncol = R, byrow = T); 
 covj = rowSums((yij - yj.mat)*(stat.mat - m.stat))/R;
 sdb = sum(covj**2)**0.5;
 return(list(m.stat = m.stat, sdb = sdb));
}

OAL.lm.smoothedboot = function(Y, A, X, R, matching.w = FALSE)
{
 ds = data.frame(Y, A, X);
 var.list = paste("X", 1:ncol(X), sep = "");
 colnames(ds) = c("Y", "A", var.list);
 bootf = function(ds, i){OAL.lm(Y = ds$Y[i], A = ds$A[i], X = ds[i,-c(1,2)], matching.w)$ATE};
 results.boot = smoothed.boot(ds, bootf, R);
 return(results.boot);
}


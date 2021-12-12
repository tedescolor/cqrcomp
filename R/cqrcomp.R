norm2 = function(v){return(sqrt(sum(v^2)))}
normInf =  function(v){return(max(abs(v)))}

Tstat = function(beta, #matrix: beta matrix values, every column for a different tau value
                 taus, # points on which the Riemann–Stieltjes is computed. The length has to be the number of rows of beta1 plus one.
                 selectedNorm = c('norm2', 'normInf')){
  beta = data.matrix(beta)
  if(ncol(beta) == 1 & length(taus) > 2 ){ # in case only 1 componet, transponse is required
    beta = t(beta)
  }
  if( ncol(beta) != length(taus) - 1){   # required for Riemann–Stieltjes integral in taus
    stop("Incorrect dimension of the inputs.")
  }
  if(selectedNorm == 'norm2'){
    currentNorm = norm2
  }else if(selectedNorm == 'normInf'){
    currentNorm = normInf
  }else{
    stop("Incorrect type.")
  }
  return( sum( apply(beta, 2, currentNorm) * ( head(taus,-1) -  tail(taus,-1)) ) )  #return  Riemann–Stieltjes of the norm of the difference of beta1 and beta2 in taus
  
}


#boot.crq function in quantreg package, but the sample object is now passed as argument instead of being computed suring the function call. 
# Therefore, the same reboostrap indeces sample can be used multiple times (in particular, for each of the two datatsets being compared.)
comp.boot.crq = function (x, y, c, taus, method, ctype = "right", R = 1, mboot, 
                          bmethod = "jack", 
                          # boostrap sample that ha sto be used, instead of be computed during the function call. 
                          # of the form sample(1:n, mboot) in case of bmethod = "jack" 
                          # of the form sample(1:n, mboot, replace = TRUE) in case bmethod = "xy-pair"
                          # of the form rexp(n) in case bmethod = "Bose"
                          bsample,
                          ...) 
{
  n <- length(y)
  p <- ncol(x)
  A <- array(0, dim = c(p, length(taus), R))
  for (i in 1:R) {
    if (bmethod == "jack") {
      s <- bsample #  bsample instead of sample(1:n, mboot)
      yb <- y[-s]
      xb <- x[-s, ]
      cb <- c[-s]
      w <- rep(1, n - mboot)
    }
    else if (bmethod == "xy-pair") {
      w <- table(bsample) # bsample instead of sample(1:n, mboot, replace = TRUE)
      s <- as.numeric(names(w))
      w <- as.numeric(w)
      yb <- y[s]
      xb <- x[s, ]
      cb <- c[s]
    }
    else if (bmethod == "Bose") {
      w <- bsample
      yb <- y
      xb <- x
      cb <- c
    }
    else stop("invalid bmethod for boot.crq")
    if (method == "Portnoy") 
      a <- crq.fit.por(xb, yb, cb, weights = w, ctype = ctype, 
                       ...)
    else if (method == "PengHuang") 
      a <- crq.fit.pen(xb, yb, cb, weights = w, ctype = ctype,...)
    else stop("Invalid method for boot.crq")
    if ((i%%floor(R/10)) == 0 & n > 1e+05) 
      cat(paste("bootstrap roughly ", 100 * (i/R), " percent complete\n"))
    A[, , i] <- coef(a, taus)
  }
  list(A = A, n = length(y), mboot = mboot, bmethod = bmethod)
}




#' Censored Quantile Regression comparison
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#' 
#' This function loads two datasets and compare the distribution of the (conditonal) survival time
#' testing the equality of the quantile regression coeffcients.
#' The datasets can be paired or independent. In the first case, they must have the same sample sizes, 
#' and the i-row of the first dataset is supposed to be correlated with i-row of the second.
#' Different statistics are possible for test: integral of the euclidean norm of the different, of the sup norm, or a Bonferroni approach testing each component separatly.
#' Differnt bootstrap method are available for the confidence interval of the test statistics. See quantreg for the details.
#' 
#' 
#' @param df1 first dataset as dataframe
#' @param df2 second dataset as dataframe
#' @param cov.names covariate names in case of conditional comparison
#' @param delta.name for delta variable (0 censored, 1 uncensored)
#' @param time.name name of the time variable
#' @param taus values of the quantiles of interest, over which compare the quantile curves. If NA, the maximum interval will be used with step 0.05, starting from 0.1
#' @param R number of bootstrap resample
#' @param verbose print also partial results,
#' @param method same as crq method, see quantreg
#' @param paired boolean indicating whether the samples are paried
#' @param test.type type of test among using norm2, normInf or bonferroni approach
#' 
#' @return pcalue of the test H0: beta1(tau) = beta2(tau), for each tau in taus, using the statistics selected in test.type, and the bootstrap confidence interval specified in bmethod

#' @export
cqrcomp = function(df1, #first dataset as dataframe
                   df2, #second dataset as dataframe
                   cov.names = NULL, #covariate names in case of conditional comparison
                   delta.name = "delta", #name for delta variable (0 censored, 1 uncensored)
                   time.name = "y", #name of the time variable
                   taus = NA, # values of the quantiles of interest, over which compare the quantile curves. 
                   R=100, #number of boostrap resample
                   verbose = F, #print also partial results,
                   method, # crq method, see quantreg
                   bmethod =c("Bose", "jack", "xy-pair"), # boostrap method
                   paired=F, #whether the samples are paried
                   test.type = c('norm2', 'normInf','bonf') # type of test among using norm2, normInf or bonferroni correction
){
  if(is.na(taus)){
    max.quantile = max(
      sum(df1[,delta.name])/nrow(df1),
      sum(df2[,delta.name])/nrow(df2)
      )
    if(max.quantile<0.1){
      stop("Censor to high, impossible to continue. Try to specify taus.")
    }
    taus = seq(0.1, max.quantile, by = 0.05)
    print("taus = "); print(taus);
    
  }
  
  if( !(bmethod %in% c("Bose", "jack", "xy-pair") ) ){
    stop("No valid boostrap method.")
  }
  n1 = nrow(df1);  n2 = nrow(df2); # save the dimension of the datasets. in case different, the asymptotic considered statistic is min(n1,n2) * (beta1-beta2)(tau)
  if(paired & n1!=n2){
    stop("Paired samples with different sizes")
  }
  if( "bonf" %in% test.type & is.null(cov.names)){
    stop("Use norm2 and not bonf if cov.names = NULL")
  }
  measurevar <- "survival::Surv(get(time.name),get(delta.name))"
  if(is.null(cov.names)){
    formula = as.formula(paste(measurevar, "1", sep=" ~ "))
  }else{ #construct formula for crq function
    groupvars  <- cov.names
    paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ ")
    formula = as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  }
  
  
  rq1 <- quantreg::crq( formula = formula ,data = df1, taus = taus, method = method)
  rq2 <- quantreg::crq( formula = formula ,data = df2, taus = taus, method = method)
  #extract the coefficients of the quantile regression estimations
  t1 = coef(rq1,taus= taus)
  t2 = coef(rq2,taus= taus)
  if(verbose){
    print(t1); print(t2);
  }
  
  if(any(is.na(t1- t2))) {
    stop("NA in the estimation. Consider different taus.")
  }
  
  #######################
  ### bootstrap procedure
  #######################
  # convert  the data in proper matrices for boostrap computation
  num.var = ifelse(is.null(cov.names),1,length(cov.names) +1)
  x1 = matrix(1, n1, num.var)
  x2 = matrix(1, n2, num.var) 
  
  if(!is.null(cov.names)){
    x1[,2:num.var] = data.matrix(df1[,cov.names])
    
    x2[,2:num.var] = data.matrix(df2[,cov.names])
    colnames(x1) <- c("intercet", cov.names)
    colnames(x2) <- c("intercet", cov.names)
  }else{
    colnames(x1) <- "intercept"
    colnames(x2) <- "intercept"
  }
  
  y1 = as.matrix(df1[,time.name])
  y2 = as.matrix(df2[,time.name])
  
  c1 = as.matrix(df1[,delta.name])
  c2 = as.matrix(df2[,delta.name])
  
  #compute bootstrap
  A1 = array(0, dim = c(ncol(x1), length(taus), 0)) #initialize bootstrap result matrix for df1 
  A2 = array(0, dim = c(ncol(x1), length(taus), 0)) #initialize bootstrap result matrix for df2
  
  #depending on the bootstrap method, different mboot. See original boot.crq for details.
  n = min(n1,n2);
  if(bmethod == "jack"){ 
    mboot <- 2 * ceiling(sqrt(n))
  }else{ mboot <- n}
  
  for(i in 1:R){
    if(paired){
      if(bmethod == "jack"){
        bsample1 = sample(1:n, mboot)
        bsample2 = bsample1
      }else if(bmethod == "xy-pair"){
        bsample1 = sample(1:n, mboot, replace = TRUE)
        bsample2 = bsample1
      }else{ #Bose method
        bsample1 = rexp(n)
        bsample2 = bsample1
      }
    }else{
      if(bmethod == "jack"){
        bsample1 = sample(1:n1, mboot)
        bsample2 = sample(1:n2, mboot)
      }else if(bmethod == "xy-pair"){
        bsample1 = sample(1:n1, mboot, replace = TRUE)
        bsample2 = sample(1:n2, mboot, replace = TRUE)
      }else{ #Bose method
        bsample1 = rexp(n1)
        bsample2 = rexp(n2)
      }
    }
    
    A1 = abind(A1, comp.boot.crq(x = x1,y = y1,c = c1,taus = taus,method = method,mboot = mboot,bmethod = bmethod,bsample = bsample1)$A,
               rev.along=1)
    A2 =abind(A2, comp.boot.crq(x = x2,y = y2,c = c2,taus = taus,method = method,mboot = mboot,bmethod = bmethod,bsample = bsample2)$A,
              rev.along = 1)
  }
  
  # constrcut Riemann–Stieltjes taus 
  if(length(taus) == 1){
    TstatTaus = c(taus[1],taus[1] + 1)
  }else{
    TstatTaus = c(taus, taus[length(taus)] + taus[length(taus)] -  taus[length(taus)-1]  )
  }
  
  result = list()
  if( "bonf" %in% test.type ){
    
    pvalues = c()
    for(i in 1:(length(cov.names)+1)){
      if(length(taus) == 1){
        pvalue =  1 - sum( apply(
          array( (A1-A2-(t1-t2))[i,,], dim = c(1,dim(A1)[2:3])), 3,
          Tstat, taus = TstatTaus,selectedNorm = "norm2" ) < Tstat( (t1-t2)[i],taus = TstatTaus, selectedNorm = "norm2" ) ) / R
        pvalues = c(pvalues, pvalue)
      }else{
        pvalue =  1 - sum( apply(
          array( (array(A1-A2, dim = dim(A1)) - array(rep((t1-t2),R),dim = c(dim(t1),R) ))[i,,], dim = c(1,dim(A1)[2:3])), 3,
          Tstat, taus = TstatTaus,selectedNorm = "norm2" ) < Tstat( (t1-t2)[i,],taus = TstatTaus, selectedNorm = "norm2" ) ) / R
        pvalues = c(pvalues, pvalue)
      }
      
    }
    pvalue = min(pvalues)
    result = c(result, list(bonf = pvalue))
  }
  if("norm2" %in% test.type){
    pvalue = 1 - sum( apply( (array(A1-A2, dim = dim(A1)) - array(rep((t1-t2),R),dim = c(dim(t1),R) )) , length(dim(A1)),Tstat, taus = TstatTaus,selectedNorm = "norm2" ) < Tstat(t1-t2,taus = TstatTaus, selectedNorm = "norm2" ) ) / R
    result = c(result, list(norm2 = pvalue))
  }
  if("normInf" %in% test.type){
    pvalue = 1 - sum( apply( (array(A1-A2, dim = dim(A1)) - array(rep((t1-t2),R),dim = c(dim(t1),R) )) , length(dim(A1)),Tstat, taus = TstatTaus,selectedNorm = "normInf" ) < Tstat(t1-t2,taus = TstatTaus, selectedNorm = "normInf" ) ) / R
    result = c(result, list( normInf = pvalue))
  }
  return(result)
}




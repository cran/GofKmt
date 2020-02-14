
#' Implementing Khmaladze Martingale Transformation.
#'
#' Performs goodness-of-fit test through Khmaladze matringale transformation
#'@param X  a random sample of n observations
#'@param strDist  the name of the null distribution for the hypothesis test: Normal, Cauchy, or Logistic. Other distributions such as gamma, Weibull and Frechet will be available in later versions. 
#'@param bGraph  a logical value which specifies whether or not to get the graph of the objective function of the martingale transformation.
#'@param nNum  the number of ticks on each segmented interval when drawing the graph of the objective function. The default is 10. Bigger value will result in a smoother graph.
#'@return A list of the following values: 
#'\describe{
#'\item{opt.x}{ the value of x where the optimum of the objective function - which is also the test statistic - occurs.} 
#'
#'\item{test.stat}{ the test statistic obtained through Khmaladze martingale transformation}
#'
#'\item{graph.data}{ a data frame which includes the information of the objective function.}
#'
#'\item{graph}{ a ggplot object which includes the graph of the objective function.}
#'
#'\item{intervals}{ a list of segmented intervals over which the graph of the objective function is defined.}
#'
#'\item{crit.value}{ a vector of critical values for the level of 0.01, 0.025, 0.05, and 0.10}
#'
#'\item{mu}{ the point estimate for the location parameter mu}
#'
#'\item{sigma}{ the point estimate for the scale parameter sigma}
#'}
#'
#'@examples
#'####################
#'n = 10
#'X = rnorm(n, 1,3)    # Generate a random sample of n observations from N(1,3)
#'strDist = "Normal"
#'lResult = KhmaladzeTrans(X, strDist, bGraph=TRUE, nNum=10)
#'KMT_OptimalX = lResult$opt.x
#'KMT_TestStat = lResult$test.stat
#'KMT_DM = lResult$graph.data
#'KMT_Graph = lResult$graph
#'KMT_Intervals = lResult$intervals
#'KMT_CriticalValue = lResult$crit.value
#'KMT_Muhat = lResult$mu
#'KMT_Sigmahat = lResult$sigma
#'#####################
#'
#'#####################
#'n = 10
#'X = rlogis(n, 1,2)  # Generate a random sample of n observations from the logistic distribution
#'strDist = "Logistic"
#'lResult = KhmaladzeTrans(X, strDist, bGraph=TRUE, nNum=10)
#'KMT_OptimalX = lResult$opt.x
#'KMT_TestStat = lResult$test.stat
#'KMT_DM = lResult$graph.data
#'KMT_Graph = lResult$graph
#'KMT_Intervals = lResult$intervals
#'KMT_CriticalValue = lResult$crit.value
#'KMT_Muhat = lResult$mu
#'KMT_Sigmahat = lResult$sigma
#'#####################
#'#####################
#'n = 10
#'X = rcauchy(n, 0,1)  # Generate a random sample of n observations from Cauchy distribution
#'strDist = "Cauchy"
#'lResult = KhmaladzeTrans(X, strDist, bGraph=TRUE, nNum=10)
#'KMT_OptimalX = lResult$opt.x
#'KMT_TestStat = lResult$test.stat
#'KMT_DM = lResult$graph.data
#'KMT_Graph = lResult$graph
#'KMT_Intervals = lResult$intervals
#'KMT_CriticalValue = lResult$crit.value
#'KMT_Muhat = lResult$mu
#'KMT_Sigmahat = lResult$sigma
#'#####################





#'@references
#'[1] Khmaladze, E.V., Koul, H.L. (2004). Martingale transforms goodness-of-fit tests in regression models. Ann. Statist., 32. 995-1034
#'@references
#'[2] E.V. Khmaladze, H.L. Koul (2009). Goodness-of-fit problem for errors in nonparametric regression: distribution free approach. Ann. Statist., 37(6A) 3165-3185.
#'@references
#'[3] Kim, Jiwoong (2020). Implementation of a goodness-of-fit test through Khmaladze martingale transformation.
#'@export
#'@useDynLib GofKmt


KhmaladzeTrans = function(X, strDist, bGraph=FALSE, nNum=10){
  

  fast.estimate = FALSE
  n = length(X)
  
  if(strDist == "Normal"){
    muhat = mean(X)
    sigmahat = sd(X)
    
  }else if(strDist == "Logistic"){
    
    if(fast.estimate == TRUE){
      muhat = median(X)
      sigmahat = sqrt(3)/pi*sd(X)
      
    }else{
      x01 = median(X)
      x02 = sqrt(3)/pi*sd(X)
      
      startVal = c(x01, x02)
      mle_method = optim(startVal, mleLogis(X))
      mle_par = mle_method$par
      
      muhat = mle_par[1]
      sigmahat = mle_par[2]
    }
    
    
  }else if(strDist == "Cauchy"){

    if(fast.estimate == TRUE){
      
      muhat = median(X)
      sighat = GetCauchyScale(X, muhat)
      
    }else{
      x01 = median(X)
      
      quan3 = quantile(X, 0.75)
      quan1 = quantile(X, 0.25)
      
      x02 = (quan3[[1]]-quan1[[1]])/2
      
      ####################
      x0 = c(x01, x02)
      LBV = c(-Inf, 0)
      UBV = c(Inf, Inf)
      
      MRe = solnp( pars=x0, fun=MLECauchyObjFunction(X), eqfun=NULL, eqB=NULL, ineqfun=NULL, ineqLB = NULL,
                   ineqUB = NULL, LB = LBV, UB = UBV)
      
      CauchyMle = MRe$pars
      
      muhat = CauchyMle[1]
      sigmahat = CauchyMle[2]
      
    }
      
    
  }else{
    message("Name of null distribution is not valid.")
    stop()
  }
  
  StandX = rep(0, times=n)
  
  for(i in 1:n){
    StandX[i] = (X[i]-muhat)/sigmahat
    
  }
  
  StandX = sort(StandX)
  
  d_alpha= c(2.80705, 2.49771, 2.24241, 1.959, 1.64)
  
  #data(Integration_Tables)
  
  NormalMat = Integration_Tables$normal.table
  LogisMat = Integration_Tables$logistic.table1
  ReMat = Integration_Tables$logistic.table2
  CauchyMat = Integration_Tables$cauchy.table
  
  lst = KmtMain(X, NormalMat, LogisMat, ReMat, CauchyMat, strDist, bGraph, nNum)
  
  lineVec = lst[[3]]
  gVec = lst[[4]]
  nLen = length(lineVec)
  
  DM = matrix(0, nLen, 3)
  DM = data.frame(DM)
  
  lineVec = round(lineVec, digits=3)
  gVec=round(gVec, digits=3)
  
  colnames(DM, do.NULL = TRUE)
  colnames(DM) = c("x", "y", "Interval")
  
  DM[, "x"] = lineVec
  DM[, "y"] = gVec
  
  nGroup = n+1
  
  xIntVec = rep(0, times=n)
  
  for(i in 1:nGroup){
    
    Sn = (i-1)*nNum+1
    En = i*nNum
    
    if(i==1){
      Sval = "-Inf"
    }else{
      Sval = DM[Sn, "x"]
    }
    
    if(i==nGroup){
      Eval = "Inf"
    }else{
      Eval = DM[En+1, "x"]
      xIntVec[i] = Eval
    }
    
    Str = paste("", i, ". [", Sval, ", ", Eval, ")", sep="")
    print(Str)
    DM[Sn:En, "Interval"] = Str
  }
  
  DM[, "Interval"] = as.factor(DM[, "Interval"])
  
  Intervals = levels(DM[, "Interval"])
  
  gObj = ggplot(data=DM, mapping=aes(x=DM[,1], y=DM[,2], group = DM[,3]))+
    geom_line(color="red")+
    labs(title = "Graph of the Objective Function", x="x", y="y")
  
  
  lResult = list(opt.x = lst[[1]], test.stat=lst[[2]], 
                 graph.data = DM, graph = gObj, intervals=Intervals, mu = muhat, sigma = sigmahat, 
                 crit.value = c("1%: 2.80705", "2.5%: 2.49771", "5%: 2.24241", "10%: 1.959") )
  return(lResult)
}


GetCauchyScale = function(XVec, mu){
  
  al = 0.5
  quant = quantile(XVec, al )
  Qval = quant[[1]]
  
  sig = (al-mu)/( tan( pi*( pcauchy(Qval, 0,1) - 0.5 ) ) )
  
  return( abs(sig)+0.0001)
}


MLECauchyObjFunction = function(XVec){
  nLength = length(XVec)
  
  Dual = function(x){
    tempsum = -nLength*log(pi*x[2])
    for(i in 1:nLength){
      tempsum = tempsum - log(1+( (XVec[i]-x[1] ) / x[2]  )^2)
    }
    return(-tempsum)
  }
  return(Dual)
}



mleLogis = function(X){
  
  nLeng = length(X)
  Dual = function(param){
    tempsum=0
    for (i in 1:nLeng){
      normx = (X[i]-param[1])/param[2]
      tempsum = tempsum + normx - log(param[2]) -2*log(1+exp(normx))
    }
    return(-tempsum)
  }
  return(Dual)
  
}













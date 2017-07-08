

#' Performs goodness-of-fit test through Khmaladze matringale transformation
#'@param X  - Random sample of n observations
#'@param F0 - Name of null distribution: Normal, Cauchy, Logistic, Gamma, Gumbel, Weibull, or Rayleigh 
#'@param Shape - The shape parameter used for Gamma and Weibull. The value should be even. 
#'@return TestStat - Test statistic obtained through Khmaladze martingale transformation 
#'@return CritValue - Vector of critical values for the level of 0.01, 0.025, 0.05, and 0.10 
#'@return Mu - Maximum likelihood estimator of location parameter mu 
#'@return Sigma - Maximum likelihood estimator of scale parameter sigma
#'@examples
#'####################
#'n = 10
#'Sample = rnorm(n, 1,3)    # Generate a random sample of n observations from N(1,3)
#'KMT_Result = KhmaladzeTrans(Sample, "Normal", Shape=2)
#'KMT_TestStat = KMT_Result$TestStat
#'KMT_CriticalValue = KMT_Result$CritValue
#'KMT_Muhat = KMT_Result$Mu
#'KMT_Sigmahat = KMT_Result$Sigma
#'#####################
#'
#'#####################
#'n = 10
#'Sample = rlogis(n, location=2, scale=3)  # Generate a random sample from logistic distribution
#'KMT_Result = KhmaladzeTrans(Sample, "Logistic", Shape=2)
#'KMT_TestStat = KMT_Result$TestStat
#'KMT_CriticalValue = KMT_Result$CritValue
#'KMT_Muhat = KMT_Result$Mu
#'KMT_Sigmahat = KMT_Result$Sigma
#'#####################



#'@references
#'[1] E.V. Khmaladze, H.L. Koul (2004). Martingale Transforms Goodness-of-fit tests in regression models. Ann. Statist., 32. 995-1034 
#'@references
#'[2] E.V. Khmaladze, H.L. Koul (2009). Goodness-of-fit problem for errors in nonparametric regression: Distribution free approach. Ann. Statist., 37(6A) 3165-3185.
#'#'@references
#'[2] H.L. Koul, L. Sakhanenko (2005). Goodness of fit testing in regression: A finite sample comparison of bootstrap methodology and Khmaladze transformation. Stat. Probab. Lett, 74 290-302. 
#'@importFrom "stats" "optim" "dnorm" "integrate" "median" "optimize" "pnorm" "qgamma" "qnorm" "quantile" "sd" 
#'@importFrom "Rsolnp" "solnp" 
#'@export


KhmaladzeTrans = function(X, F0, Shape=2){
  
  k=Shape
  
  fl = floor(k/2)
  resid = k-fl*2
  
  if((k<= 0)){
    message("Shape should be positive.")
    stop()
  }else if( (resid != 0) ){
    message("Shape should be even.")
    stop()
  }
  
  SampleSize = length(X)
  
  if(F0 == "Normal"){
    muhat = mean(X)
    sigmahat = sd(X)
    
  }else if(F0 == "Cauchy"){
    
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
    
  }else if(F0 == "Logistic"){
    
    x01 = median(X)
    x02 = sqrt(3)/pi*sd(X)
    
    startVal = c(x01, x02)
    mle_method = optim(startVal, mleLogis(X))
    mle_par = mle_method$par
    
    muhat = mle_par[1]
    sigmahat = mle_par[2]
    
  }else if(F0 == "Gumbel"){
    EMC = 0.5772
    sigmahat = sqrt(6/pi^2) *sd(X)
    muhat = mean(X)-sigmahat*EMC
    
  }else if(F0 == "Weibull"){
    
    muhat=0
    sigmahat = mean(X)/gamma(1+1/k)
    
  }else if(F0 == "Gamma"){
    
    muhat=0
    sigmahat = mean(X)/k
    
  }else if(F0 == "Rayleigh"){
    
    muhat=0
    sigmahat = mean(X)/sqrt(pi/2)
    
  }else{
    message("Name of null distribution is not valid.")
    stop()
  }
  
  StandX = rep(0, times=SampleSize) 
  
  for(i in 1:SampleSize){
    StandX[i] = (X[i]-muhat)/sigmahat
    
  }
  
  d_alpha= c(2.80705, 2.49771, 2.24241, 1.959, 1.64)
  
  if(F0 == "Normal"){
    KhmResult = SupWNormal(StandX)
  }else if(F0 == "Cauchy"){
    KhmResult = SupWCauchy(StandX)
  }else if(F0 == "Logistic"){
    KhmResult = SupWLogistic(StandX)
  }else if(F0 == "Gumbel"){
    KhmResult = SupWGumbel(StandX)
  }else if(F0 == "Weibull"){
    KhmResult = SupWWeibull(StandX, k)
  }else if(F0 == "Gamma"){
    KhmResult = SupWGamma(StandX, k)
  }else if(F0 == "Rayleigh"){
    KhmResult = SupWRayleigh(StandX)
  }
  
  
  lst = list(TestStat=KhmResult$objective, CritValue = c("1%: 2.80705", "2.5%: 2.49771", "5%: 2.24241", "10%: 1.959"), 
             Mu = muhat, Sigma = sigmahat )
  return(lst)
}


########################### Rayleigh Part

##############distribution

RaF = function(x){
  
  val = 1 - exp(-x^2/2)
  return(val)
  
}


###### density
Raf = function(x){
  val = x*exp(-x^2/2)
  return(val)
}

######## derivative of density

Raf1 = function(x){
  val = exp(-x^2/2) * (1-x^2)
  return(val)
}

RaIntGam22 = function(x){
  
  val = ( 1/x -2*x + x^3 )*exp(-x^2/2)
  return(val)
}


RaGam22 = function(x){
  
  IntResult = integrate(RaIntGam22, x, Inf)
  val = IntResult$value
  return(val)
}



################# det of Info Matrix

Radet = function(x){
  
  Gam11 = 1 - RaF(x)
  Gam12 = -Raf(x)
  Gam22 = RaGam22(x)
  
  val = Gam11*Gam22 - (Gam12)^2
  return(val)
  
}


################# 1/det of A22

Ra1det = function(x){

  Gam11 = 1 - RaF(x)
  Gam12 = -Raf(x)
  Gam22 = RaGam22(x)
  
  TempVal = Gam11*Gam22 - (Gam12)^2  
  val = 1/TempVal
  return(val)
  
}



######################### l(Xi)'Gam(y)l(y)f(y)

ItgofIntRayleigh = function(Xi){
  
  dual = function(y){
    
    det = Radet(y)
    Gam22 = RaGam22(y)
    
    Fy = RaF(y)
    f = Raf(y)
    f1 = Raf1(y)
    f2 = f^2 
    fxi = Raf(Xi)
    f1xi = Raf1(Xi)
    
    val = ( fxi*f*Gam22 + fxi*f*f1 + f1xi*f2+f1xi*f1-f1xi*f1*Fy )/(fxi*det)
    
    return(val)
    
  }
  return(dual)
}



IntRayleigh = function(Xi,t){
  Ft1 = sqrt( -2*log(1-t) ) 
  upValVec = c(Xi, Ft1)
  upVal = min(upValVec)
  ans = integrate( ItgofIntRayleigh(Xi), lower=0, upper = upVal, subdivisions = 2000)$value
  return(ans)
  
}


WRayleigh = function(X){
  
  InsideW1n = function(t){
    z = sqrt( -2*log(1-t) )
    n=length(X)
    tempsum = 0
    for (i in 1:n){
      
      Xi=X[i]
      if(Xi <= z){Ixi=1}
      else{Ixi=0}
      tempsum = tempsum + ( Ixi - IntRayleigh(Xi,t)  )
      
    }
    tempsum = tempsum/sqrt(n)
    
    AbsVal = abs(tempsum)
    
    if(AbsVal> 1000){return(1000)}
    
    return(AbsVal)	
  }
  return(InsideW1n)
  
}

SupWRayleigh = function(X){
  
  tmin = optimize(WRayleigh(X),  lower=0, upper =1, maximum=TRUE )
  
}



#####################################




####################### Gamma Part

##############distribution
IntofGF = function(k){
  
  dual = function(x){
    val = x^(k-1)*exp(-x)  
    return(val)
  }
  return(dual)
}

GF = function(k){
  
  dual = function(x){
    val = integrate(IntofGF(k), 0, x)$value
    return(val/gamma(k))
  }
  return(dual)
}


###### density
Gf = function(k){
  dual = function(x){
    val = x^(k-1)*exp(-x)/gamma(k)
    return(val)
  }
  return(dual)
}

######## derivative of density

Gf1 = function(k){
  dual = function(x){
    val = x^(k-2)*exp(-x)*( (k-1)-x)/gamma(k)
    return(val)
  }
  return(dual)
} 

GIntGam22 = function(k){
  dual = function(x){
    val = ( (k-1)^2/x^2 - 2*(k-1)/x + 1 )*x^(k-1)*exp(-x)/gamma(k)
    return(val)
  }
  return(dual)
}

GGam22 = function(k){
  
  dual = function(x){
    
    IntResult = integrate(GIntGam22(k), x, Inf)
    val = IntResult$value
    return(val)
  }
  return(dual)  
  
}



################# det of Info Matrix

Gdet = function(k){
  
  dual = function(x){
    
    Gam11 = 1 - GF(k)(x)
    Gam12 = -Gf(k)(x)
    Gam22 = GGam22(k)(x)
    
    val = Gam11*Gam22 - (Gam12)^2
    return(val)
  }
  return(dual)
}


################# 1/det of A22

G1det = function(k){
  
  dual = function(x){
    
    Gam11 = 1 - GF(k)(x)
    Gam12 = -Gf(k)(x)
    Gam22 = GGam22(k)(x)
    
    TempVal = Gam11*Gam22 - (Gam12)^2
    val = 1/TempVal
    return(val)
  }
  return(dual)
}



######################### l'Gam(y)l(y)f(y)


ItgofIntGamma = function(Xi, k){
  
  dual = function(y){
    
    det = Gdet(k)(y)
    Gam22 = GGam22(k)(y)
    
    Fy = GF(k)(y)
    f = Gf(k)(y)
    f1 = Gf1(k)(y)
    f2 = f^2 
    fxi = Gf(k)(Xi)
    f1xi = Gf1(k)(Xi)
    
    val = ( fxi*f*Gam22 + fxi*f*f1 + f1xi*f2+f1xi*f1-f1xi*f1*Fy )/(fxi*det)
    
  }
  return(dual)
}



IntGamma = function(Xi,t, k){
  upValVec = c(Xi, qgamma(t, shape=k, scale = 1))
  upVal = min(upValVec)
  ans = integrate( ItgofIntGamma(Xi, k), lower=0, upper = upVal, subdivisions = 3000)$value
  return(ans)
  
}


WGamma = function(X, k){
  
  InsideW1n = function(t){
    z = qgamma(t, shape=k, scale = 1)
    n=length(X)
    tempsum = 0
    for (i in 1:n){
      
      Xi=X[i]
      if(Xi <= z){Ixi=1}
      else{Ixi=0}
      tempsum = tempsum + ( Ixi - IntGamma(Xi,t, k)  )
      
    }
    tempsum = tempsum/sqrt(n)
    
    AbsVal = abs(tempsum)
    
    if(AbsVal> 1000){return(1000)}
    
    return(AbsVal)	
  }
  return(InsideW1n)
  
}

SupWGamma = function(X, k){
  
  tmin = optimize(WGamma(X, k),  lower=0, upper =1, maximum=TRUE )
  
}



#####################################






###########################
######################### Weibull Part

##############distribution
WBF = function(k){
  
  dual = function(x){
    val = 1-exp(-x^k)  
    return(val)
  }
  return(dual)
}


###### density
WBf = function(k){
  dual = function(x){
    val = k*x^(k-1)*exp(-x^k)
    return(val)
  }
  return(dual)
}

######## derivative of density

WBf1 = function(k){
  dual = function(x){
    val = k*x^(k-2)*exp(-x^k)*( (k-1)-k*x^k )
    return(val)
  }
  return(dual)
} 

WBIntGam22 = function(k){
  dual = function(x){
    val = exp(-x)*( (k-1)*x^(-1/k) - k*x^((k-1)/k)  )^2
    return(val)
  }
  return(dual)
}

WBGam22 = function(k){
  
  dual = function(x){
    xk = x^k
    IntResult = integrate(WBIntGam22(k), xk, Inf)
    val = IntResult$value
    return(val)
  }
  return(dual)  
  
}

################# det of Info Matrix

WBdet = function(k){
  
  dual = function(x){
    
    Gam11 = 1 - WBF(k)(x)
    Gam12 = -WBf(k)(x)
    Gam22 = WBGam22(k)(x)
    
    val = Gam11*Gam22 - (Gam12)^2
    return(val)
  }
  return(dual)
}


################# 1/det of Info Matrix

WB1det = function(k){
  
  dual = function(x){
    
    Gam11 = 1 - WBF(k)(x)
    Gam12 = -WBf(k)(x)
    Gam22 = WBGam22(k)(x)
    
    TempVal = Gam11*Gam22 - (Gam12)^2
    val = 1/TempVal
    return(val)
  }
  return(dual)
}


######################### l'Gam(y)l(y)f(y)


ItgofIntWeibull = function(Xi, k){
  
  dual = function(y){
    
    det = WBdet(k)(y)
    Gam22 = WBGam22(k)(y)
    
    Fy = WBF(k)(y)
    f = WBf(k)(y)
    f1 = WBf1(k)(y)
    f2 = f^2 
    fxi = WBf(k)(Xi)
    f1xi = WBf1(k)(Xi)
    
    val = ( fxi*f*Gam22 + fxi*f*f1 + f1xi*f2+f1xi*f1-f1xi*f1*Fy )/(fxi*det)
    
  }
  return(dual)
}



IntWeibull = function(Xi,t, k){

  Ft1 = (-log(1-t))^(1/k)
  upValVec = c(Xi, Ft1)
  upVal = min(upValVec)
  ans = integrate( ItgofIntWeibull(Xi, k), lower=0.01, upper = upVal, subdivisions=2000)$value
  return(ans)
  
}


WWeibull = function(X, k){
  
  InsideW1n = function(t){
    z = (-log(1-t))^(1/k)
    n=length(X)
    tempsum = 0
    for (i in 1:n){
      
      Xi=X[i]
      if(Xi <= z){Ixi=1}
      else{Ixi=0}
      tempsum = tempsum + ( Ixi - IntWeibull(Xi,t, k)  )
      
    }
    tempsum = tempsum/sqrt(n)
    
    AbsVal = abs(tempsum)
    
    if(AbsVal> 1000){return(1000)}
    
    return(AbsVal)	
  }
  return(InsideW1n)
  
}

SupWWeibull = function(X, k){
  
  tmin = optimize(WWeibull(X, k),  lower=0, upper =1, maximum=TRUE )
  
}



#####################################



####################### Gumbel Part


wy = function(y){
  exp(-y)
}

Dety = function(y){
  (1-exp(-wy(y))) * (1/(wy(y))^2 - (1/(wy(y))^2+1)*exp(-wy(y))  ) - exp(-2*wy(y))
}



ItgofIntGumbel = function(z){
  
  dual = function(y){
    Numer =  exp(-wy(y))/wy(y) - exp(-2*wy(y))/wy(y)- exp(-2*wy(y)) +
      (-1+exp(-z)) * (exp(-wy(y)) - exp(-wy(y))/wy(y) + exp(-2*wy(y))/wy(y)  )
    val = Numer/Dety(y)
    return(val)
    
  }
  return(dual)
}



IntGumbel = function(ri,t){
  Ft1 = -log(-log(t))
  upValVec = c(ri, Ft1)
  upVal = min(upValVec)
  ans = integrate( ItgofIntGumbel(ri), lower=-300, upper = upVal)$value
  return(ans)
  
}


WGumbel = function(r){
  
  InsideW1n = function(t){
    x = -log(-log(t))
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      
      tempsum = tempsum + ( (r[i]<=x) - IntGumbel(r[i],t)  )
      
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWGumbel = function(r){
  
  tmin = optimize(WGumbel(r),  lower=0, upper =1, maximum=TRUE )
  
}



#############################


################### Normal Part
afy = function(y){
  dnorm(y, 0,1)/(1-pnorm(y,0,1))
} 

ItgofIntNormal = function(Xi){
  
  dual = function(y){
    a = afy(y)
    
    val = ( a - Xi*(a-y)*a )/ ( y*a+1-a^2) 
    return(val)
  }
  return(dual)
}


IntNormal = function(ri,t){
  upValVec = c(ri, qnorm(t, 0,1))
  upVal = min(upValVec)
  ans = integrate( ItgofIntNormal(ri), lower=-Inf, upper = upVal)$value
  return(ans)
  
}


WNormal = function(r){
  
  InsideW1n = function(t){
    x = qnorm(t, 0,1)
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      
      tempsum = tempsum + ( (r[i]<=x) - IntNormal(r[i],t)  )
      
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWNormal = function(r){
  
  tmin = optimize(WNormal(r),  lower=0, upper =1, maximum=TRUE )
  
}


##############################################



#################################### Cauchy

##################### Cauchy MLE

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


#####################################

################################# Cauchy

##############distribution
CauF = function(x){
  
  val = 0.5+atan(x)/pi  
  return(val)
}


###### density
Cauf = function(x){
  val = (1+x^2)^(-1)/pi
  return(val)
}

######## derivative of density

Cauf1 = function(x){
  val = (-2*x)*(1+x^2)^(-2)/pi
  return(val)
}

CauGam22 = function(x){
  
  val = pi/4+ x/(1+x^2)^2-0.5*atan(x) - 0.5*x/(1+x^2)
  
  val = val/pi
  return(val)
}



################# det of Info Matrix

Caudet = function(x){
  
  Gam11 = 1 - CauF(x)
  Gam12 = -Cauf(x)
  Gam22 = CauGam22(x)
  
  val = Gam11*Gam22 - (Gam12)^2
  return(val)
  
}


################# 1/det of A22

Cau1det = function(x){
  
  Gam11 = 1 - CauF(x)
  Gam12 = -Cauf(x)
  Gam22 = CauGam22(x)
  
  TempVal = Gam11*Gam22 - (Gam12)^2
  
  val = 1/TempVal
  return(val)
  
}



######################### l(Xi)'Gam(y)l(y)f(y)

ItgofIntCauchy = function(Xi){
  
  dual = function(y){
    
    det = Caudet(y)
    Gam22 = CauGam22(y)
    
    Fy = CauF(y)
    f = Cauf(y)
    f1 = Cauf1(y)
    f2 = f^2 
    fxi = Cauf(Xi)
    f1xi = Cauf1(Xi)
    
    val = ( fxi*f*Gam22 + fxi*f*f1 + f1xi*f2+f1xi*f1-f1xi*f1*Fy )/(fxi*det)
    
    return(val)
    
  }
  return(dual)
}


KhmIntValCauchy = function(ri,t){
  
  Ft1 = tan(pi*(t-1/2))
  upValVec = c(ri, Ft1)
  
  upVal = min(upValVec)
  
  funcVal = integrate( ItgofIntCauchy(ri), lower=-100, upper = upVal)$value
  return(funcVal)
}


###################  


WCauchy = function(r){
  
  InsideW1n = function(t){
    x = tan(pi*(t-1/2))
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      
      tempsum = tempsum + ( (r[i]<=x) - KhmIntValCauchy(r[i],t)  )
      
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWCauchy = function(r){
  
  tmin = optimize(WCauchy(r),  lower=0, upper =1, maximum=TRUE )
  
}



########################### Logistic

#######################

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
###########

KhmIntValLogistic = function(ri,t){
  
  Ft1 = -log(1/t-1)
  upValVec = c(ri, Ft1)
  upVal = min(upValVec)
  
  val = 6*exp(upVal)/(1+exp(ri)) - 2*log( 1+ exp(upVal)) 
  return(val)
}

######################

WLogistic = function(r){
  
  InsideW1n = function(t){
    x = -log(1/t-1)
    n=length(r)
    tempsum = 0
    for (i in 1:n){
      tempsum = tempsum + ( (r[i]<=x) - KhmIntValLogistic(r[i],t)  )
    }
    tempsum = tempsum/sqrt(n)
    return(abs(tempsum))	
  }
  return(InsideW1n)
  
}

SupWLogistic = function(r){
  
  tmin = optimize(WLogistic(r),  lower=0, upper=1, maximum=TRUE )
  
}









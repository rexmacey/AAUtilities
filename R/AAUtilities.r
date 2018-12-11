#' Diversification score
#' 
#' The diversification score for a portfolio of n equally-weighted assets will
#' be n.  All things being equal, the larger the value the greater the
#' diversification. This is a simple measure ignoring the covariance matrix. Use
#' with caution.
#' 
#' @param w Weights of a portfolio
#'   
#' @return value Diversification score
#' @export
#' 
#' @examples diversification_score(c(.5,.4,.1))
diversification_score<-function(w){
    out<-1/sum(w^2)
    return(out)
}

ConfidenceIntervalAnnualized<-function(p,ExpRetAnnual,SD,T,method="RAFI"){
    z<-qnorm(p,0,1)
    out<-ExpRetAnnual+z*SD/(2*sqrt(T)) #per page 16 of RAFI AA-Asset_Class-Risk.pdf
    return(out)
}

#' Future Wealth
#' 
#' Calculates the wealth for which there is only a p chance of doing worse than
#' over t years given a mean of r and a standard deviation of sd. Assumes a
#' lognormal distribution.  Specify inputs in percentage formats (e.g. 8%=8, not
#' 0.08).
#' @param p Probability of doing worse than result
#' @param r Mean return
#' @param sd Standard deviation of returns
#' @param t Period in years
#' @param amt Iniitial investment. Default = 1
#'   
#' @return value Future wealth
#' @export
#' 
#' @examples AA_Wealth(.5,8,12,10)
#' 
AA_Wealth<-function(p,r,sd,t,amt=1){
    vmean<-1+r/100
    vsd<-sd/100
    vlnsd <- sqrt(log(1 + (vsd / vmean) ^ 2))
    vlner <- log(vmean) - vlnsd ^ 2 / 2
    return(amt * exp(qnorm(1 - p, vlner * t, vlnsd * sqrt(t))))
}

#' Lognormally distributed return
#' 
#' Calculates the return for which there is only a p chance of doing worse
#' than over t years given a mean of r and a standard deviation of
#' sd. Periods over 1 year are annualized. Specify inputs in percentage formats (e.g. 8%=8, not 0.08).
#' @param p Probability of doing worse than result
#' @param r Mean return
#' @param sd Standard deviation of returns
#' @param t Period in years
#'
#' @return value Return
#' @export
#'
#' @examples AA_Return(.5,8,12,10)
#' 
AA_Return<-function(p,r,sd,t){
    x <- AA_Wealth(p, r, sd, t)
    if (t<=1){
        out<-(x-1)*100
    } else {
        out<-(x^(1/t)-1)*100
    }
    return(out)
}

#' Probability of achieving a specified return
#' 
#' Returns the probability of doing at least as well as  over t
#' 
#' @param target Target (specified) return
#' @param r Mean return
#' @param sd Standard deviation of returns
#' @param t Time in years
#'   
#' @return value representing a probability.
#' @export
#' 
#' @examples AA_Prob(7.339443,8,12,10)
#' 
AA_Prob<-function(target,r,sd,t){
    vmean<-1+r/100
    vsd<-sd/100
    vlnsd <- sqrt(log(1 + (vsd / vmean) ^ 2))
    vlner <- log(vmean) - vlnsd ^ 2 / 2
    return(1 - pnorm(log((target / 100 + 1) ^ t), vlner * t, vlnsd * sqrt(t)))
}

#' Random return(s) from a lognormally distribution
#' 
#' Returns n lognormally distributed random variables with a mean of r and a
#' standard deviation of sd.  Periods over 1 year are annualized.  Specify
#' inputs in percentage formats (e.g. 8%=8, not 0.08).
#' 
#' @param n Number of returns to generate
#' @param r Mean return
#' @param sd Standard deviation of returns
#' @param t Time in years
#' @param seed Random seed
#'   
#' @return value Random returns
#' @export
#' 
#' @examples AA_Rand(10,8,12,10)
#' 
AA_Rand<-function(n,r,sd,t,seed=NA){
    if (!is.na(seed)) {
        set.seed(seed)    
    }
    vMean <- 1 + r / 100
    vSD <- sd / 100
    vLNSD <- sqrt(log(1 + (vSD / vMean) ^ 2)) # var
    vLNER <- log(vMean) - vLNSD ^ 2 / 2
    dblRnd<-runif(n)
    out <- exp(qnorm(dblRnd, vLNER * t, vLNSD * sqrt(t)) / t)
    #out<-100 * (exp(qnorm(dblRnd, vLNER * t, vLNSD * sqrt(t)) / t) - 1)
    if(t < 1){
      out <- 100 * (out^t - 1)
    } else {
      out <- 100 * (out - 1)
    }
    return(out)
}

cf.npv<-function(x,cf){
    n<-length(cf)
    d<-rep(1/(1+x),n)
    d<-d^(seq(0,n-1))
    return(sum(cf*d))
}

cf.irr<-function(cf){
    ur<-uniroot(cf.npv,c(-.99,.99),cf=cf)
    return(ur$root)
}

mc_trial<-function(r,cf,inflation=0){ # cf should be 1 longer than r
    # returns a vector of values applying a vector of returns to a vector of cash flows  
    out<-rep(0,length(cf))
    rp1<- 1+r
    infladj<-(1+inflation)^seq(0,length(cf)-1)
    out[1]<- -cf[1]
    for (i in 2:length(cf)){
        if (out[i-1]<=0) {
            break
        } else {
            out[i]<-out[i-1]*rp1[i-1] - cf[i]
            out[i]<-max(out[i],0)
        }
    }
    return(out/infladj)
}

mc_trials<-function(r,sd,trials,cf,inflation=0){
    n<-(length(cf)-1)*trials
    rmat<-matrix(AA_Rand(n,r*100,sd*100,1,101),nrow=length(cf)-1,ncol=trials)/100
    rownames(rmat)<-seq(1,length(cf)-1)
    out<-list()
    out$rmat<-rmat
    mcs<-apply(rmat,2,mc_trial,cf=cf,inflation=inflation)
    rownames(mcs)<-seq(0,length(cf)-1)
    out$trials<-mcs
    return(out)
}

psuccess<-function(v,goal){
    return(sum(v>=goal)/length(v))
}

mc_simulate<-function(r,sd,trials=1000,cf,inflation=0, t=c(1,5,10,20,Inf),
                      probs=c(.05,.25,.5,.75,.95),goal=.01){
    out<-mc_trials(r,sd,trials,cf,inflation)
    tidx<-unique(sapply(t+1,min,length(cf)))
    x<-out$trials[tidx,]
    out$percentiles<-apply(x,1,quantile,probs=probs)
    out$means<-apply(x,1,mean)
    out$psuccess<-apply(x,1,psuccess,goal=goal)
    out$assumed.return<-r
    out$assumed.sd<-sd
    out$ntrials=trials
    out$cf<-cf
    out$inflation<-inflation
    out$goal<-goal
    out$obs.return.mean<-mean(out$rmat)
    out$obs.return.sd<-sd(out$rmat)
    return(out)
}

Get10YrBEInflationRate<-function(dt=format(Sys.Date(),"%Y-%m-01")){
    library(Quandl)
    #Quandl.api_key(Sys.getenv("QUANDL_TOKEN"))
    x<-Quandl("FRED/T10YIE",start_date = "2015-01-01",end_date = dt) # for daily
    #x<-Quandl("FRED/T10YIEM",start_date = "2015-01-01",end_date = dt) # for monthly
    out<-list()
    out$date<-x[1,1]
    out$rate<-x[1,2]
    return(out)
}
  

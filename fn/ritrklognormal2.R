# functions for fitting the Ricker SR curve
# all code by NJB, Mar.19/93
# modified for R by CM
# reparameterized by KT

#KT
srfrk <- function(p1,p2,S) {
# Stock Recruitment Function: Ricker
#  return(alpha*S*exp(-((alpha*S)/(2.178*Rmax))))
# see Brander and Mohn 2004 and erratum
  return(2.71828*p2/p1*(S)*exp(-S/p1))
}

srfrkpv <- function(pv,S) {
# Stock Recruitment Function (parameter vector version): Ricker
  return(srfrk(pv[1],pv[2],S))
}

# initial parameter guess for Ricker
rkinit <- function(s,r) {
# Calculate initial parameters for fitting re-formulated Ricker
  middle <- length(r)/2
  orderRec <- order(r)
  sortedRec <- r[orderRec]
  recsortedStock <- s[orderRec]
  Rmax <- recsortedStock[middle]
  medRec <- sortedRec[middle]
#  a <- 2*medRec/Rmax
  a <- s[length(s)/2]
  param <- c(Rmax,medRec)
  return(param)
}

lnlnegloglik <- function(fitted,observed) {
  # NJB, July 27, 1994
  # Lognormal [lnl] negative [neg] log likelihood [lik]
  # Arguments:
  #   fitted   : fitted median recruitment (assumed to be free of missing values)
  #   observed : observed recruitment (assumed to be free of missing values)
  len <- length(fitted)
  sumlogr <- sum(log(observed))
  sumsqrlogratio <- sum((log(observed/fitted))^2)  
  result <- 0.5*len*log( (2*pi/len)*sumsqrlogratio) + sumlogr + len/2
  return(result)  
}

ml.rklnl <- function(s,r,ip,max.fcal=200,max.iter=20000) {
  choice <- !is.na(r) & !is.na(s)
  r <- r[choice]
  s <- s[choice]
  assign(".Stock",s)       # Store in expression frame
  assign(".Recruit",r)     # Store in expression frame
  if (missing(ip)) {
    ip <- rkinit(.Stock,.Recruit)
  }

  logalpha <- log(ip[1])
  logbeta  <- log(ip[2])
  logip <- c(logalpha,logbeta)

  rk.nll <- function(x) {
    fitted <- srfrkpv(exp(x),.Stock)
    return(lnlnegloglik(fitted,.Recruit))  
  }          
             
  nlmin.out <- optim(par=logip, fn=rk.nll,control=list(maxit=max.iter),method="Nelder-Mead", hessian=TRUE)
#  nlmin.out <- optim(par=logip, fn=rk.nll,control=list(maxit=max.iter),method="SANN", hessian=TRUE)
  se.x<-sqrt(diag(solve(nlmin.out$hessian)))
  x <- exp(nlmin.out$par)
  up.x<-exp(nlmin.out$par+1.96*se.x)
  low.x<-exp(nlmin.out$par-1.96*se.x)  
  fitted <- srfrkpv(x,.Stock)
  nll <- lnlnegloglik(fitted,.Recruit)
  sigma.squared <- sum( (log(.Recruit/fitted))^2 )/length(.Recruit)
  x[1] <- x[1]*exp(sigma.squared/2)
  up.x[1]<-up.x[1]*exp(sigma.squared/2)
  low.x[1]<-low.x[1]*exp(sigma.squared/2)
  return(list(pv=x,upper=up.x, lower=low.x, sigma=sqrt(sigma.squared),nll=nll, converged=nlmin.out$convergence,conv.message=nlmin.out$message))
}

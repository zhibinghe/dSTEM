#' Generate simulated signals
#'
#' @param l length of data, if data is periodic then the length in each period
#' @param h numerical vector of true change point locations
#' @param jump numerical vector of jump size at change point locations
#' @param b1 numerical vector of piecewise slopes
#' @param rep number of periods if data is periodic, default is 1
#' @param shift numerical vector of vertical shifts for each period, default is 0
#'
#' @return a vector of simulated signal
#' @export
#'
#' @examples
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = rep(0,7)
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5)/50
#' beta1 = c(beta1,-sum(beta1*(c(h[1],diff(h))))/(l-tail(h,1)))
#' signal = gen.signal(l,h,jump,beta1)
gen.signal = function(l, h, jump, b1, rep=1, shift=0){
  if (length(rep)!=length(shift)) stop("rep and shift should have the same length")
  if (h[length(h)] > l) stop("h should not be greater than l")
  if (length(jump)!=length(b1) | length(jump)!= length(h)+1)
    stop("Length of jump or b1 is not matched with the length of h")
  f = function(b0){
    s = vector();t = 1:l
    for(i in 1:(length(h)+1)){
      nn = ifelse(i==1,h[i],ifelse(i==(length(h)+1),l-h[i-1],h[i]-h[i-1]))
      start = ifelse(i==1,1,h[i-1]+1)
      end = ifelse(i==(length(h)+1),l,h[i])
      temp = rep(b0[i],nn) + rep(b1[i],nn)*t[start:end]
      s = append(s,temp)
    }
    return(s)
  }
  b0 = cumsum(c(jump[1],-diff(b1)*h+jump[-1]))
  return(rep(f(b0),rep)+rep(shift,each=l))
}

#' Smoothing data using Gaussian kernel
#'
#' @param x numeric vector of values to smooth
#' @param gamma bandwidth of Gaussian kernel
#'
#' @return vector of smoothed values
#' @export
#' @examples
#' smth.gau(x=rnorm(1000), gamma=20)
smth.gau = function(x, gamma){
  # Gaussian kernel
  .kern = function(x,v=1){
    temp = exp(-x^2/(2*v^2))
    return(temp/sum(temp))}
  k = ifelse(2*6*gamma <= 0.9*length(x),6,floor(0.9*length(x)/(2*gamma))) #6sigma
  Lwindow = (-k*gamma):(k*gamma)
  w = .kern(Lwindow,v=gamma)
  sx = conv(x,w,"same")
  # adjusted weights
  adj.w = function(w){
    hw = floor(length(w)/2)
    a = cumsum(w)[(hw+1):length(w)]
    return(c(a,rep(1,length(x)-2*length(a)),rev(a)))
  }
  return(sx/adj.w(w))
}

#' Estimate variance of smoothed Gaussian noise
#' @description   Estimate variance of smoothed Gaussian noise through its second-order derivative
#' @param x numerical vector of second-order derivative of kernel smoothed data
#' @param gamma bandwidth of Gaussian kernel
#' @param k numerical value, local maxima (minima) are presumed beyond \eqn{Mean(x) ± k*SD(x)}
#'
#' @return value of estimated variance of smoothed noise
#' @export
#' @examples
#' l=15000; h = seq(150,l,150)
#' jump = rep(0,length(h)+1); b1 = seq(from=0,by=0.15,length = length(h)+1)
#' signal = gen.signal(l,h,jump,b1)
#' data = signal + rnorm(length(signal),0,1) # standard white noise
#' gamma = 10
#' ddy = diff(smth.gau(data,gamma),differences=2)
#' est.sigma2(ddy,gamma,k=0.5) # true value is \eqn{\frac{1}{2\sqrt{pi}\gamma}}
est.sigma2 = function(x, gamma, k=0.5){
  lmax = which.peaks(x)
  lmin = which.peaks(x, decreasing=T)
  J = c(x[lmax],-x[lmin])
  ind = c(lmax,lmin)[J>= mean(J)+k*sd(J)]
  delete = function(y){
    v=vector()
    for(i in 1:length(y)){
      v1 = (y[i]-4*gamma) : (y[i]+4*gamma)
      v = append(v,v1)
    }
    intersect(1:length(x),unique(v))
  }
  var2d = var(x[-delete(ind)])
  return(1/(2*sqrt(pi)*(3/(8*sqrt(pi)*var2d))^(1/5)))
}

#' Estimate piecewise slope for piecewise linear model
#'
#' @param x numerical vector of signal-plus-noise data
#' @param breaks numerical vector of change-point locations
#'
#' @return a vector of estimated piecewise slope
#' @import MASS
#' @export
est.slope = function(x, breaks){
  if(length(breaks)==0){
    t = 1:length(x)
    slope = coef(lm(x~t))[2]
  }
  else{
    breaks = sort(breaks)
    k=length(breaks)
    slope = vector(length=k+1)
    for(i in 1:(k+1)){
      if(i==1) seg = 1:breaks[i]
      else if(i==k+1) seg = (breaks[i-1]+1):length(x)
      else {seg=(breaks[i-1]+1):breaks[i]}
      slope[i] = as.numeric(coef(MASS::rlm(y~t,data.frame(t=seg,y=x[seg])))[2])
    }
  }
  return(slope)
}

#' Identify pairwise local maxima and local minima of the second-order derivative
#'
#' @param vall vector of locations of significant local minima
#' @param peak vector of locations of significant local maxima
#' @param gamma bandwidth of Gaussian kernel smoothing function
#'
#' @return a list of detected pairs and detected change-point locations through second-order derivative testing
#' @export
#'
est.pair = function(vall, peak, gamma){
  # if paired, |peak-vall| should be around 2gamma
  fun = function(x,y){
    l1 = x-2.5*gamma
    l2 = x-1.5*gamma
    r1 = x+1.5*gamma
    r2 = x+2.5*gamma
    t = y[which((y>=l1&y<=l2) | (y>=r1&y<=r2))]
    if(length(t)==0) t = NA
    if(length(t)>1) t = t[which.min(abs(t-x))]
    return(t)
  }
  if(length(vall)==0) pair = NULL
  else pair = sapply(vall,fun,y=peak,simplify=TRUE)
  single = setdiff(peak,pair)
  pairs = as.data.frame(rbind(c(vall,single),c(pair,rep(NA,length(single)))))
  return(list(pair=pairs,cp=sort(as.integer(colMeans(pairs,na.rm=TRUE)))))
}

#' Multiple testing of change points for kernel smoothed data
#'
#' @param x vector of kernel smoothed data
#' @param order order of derivative of data
#' @param alpha global significant level
#' @inheritParams smth.gau
#' @param sigma standard deviation of kernel smoothed noise
#' @param breaks vector of rough estimate of change-point locations, only required when order is 1.
#' @param slope vector of rough estimate of slopes associated with \code{breaks}, only required when order is 1.
#' @param untest vector of locations unnecessary to test
#' @param nu standard deviation of Gaussian kernel used to generate autocorrelated Gaussian noise,
#' it is 0 if the noise is Gaussian white noise.
#' @param is.constant logical value indicating if the signal is piecewise constant,
#' if TRUE, \code{breaks} and \code{slope} are not necessary.
#' @param margin length of one period of data \code{x}
#'
#' @return a list of estimated change-point locations and threshold for p-value
#' @export
#' @examples
#' ## piecewise linear signal
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = rep(0,7)
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5)/50
#' beta1 = c(beta1,-sum(beta1*(c(h[1],diff(h))))/(l-tail(h,1)))
#' signal = gen.signal(l,h,jump,beta1)
#' noise = rnorm(length(signal),0,2)
#' gamma = 25
#' sdata = smth.gau(signal+noise,gamma)
#' ddy = diff(sdata,differences=2)
#' model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05)
#' ## piecewise constant
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)
#' beta1 = rep(0,length(h)+1)
#' signal = gen.signal(l,h,jump,beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' sdata = smth.gau(signal+noise,gamma)
#' dy = diff(sdata)
#' model1 = cpTest(x=dy,order=1,alpha=0.05,gamma=gamma,is.constant=TRUE)
#' ## piecewise linear with jump
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)*3
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5,-0.5)/50
#' signal = gen.signal(l=l,h=h,jump=jump,b1=beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' sdata = smth.gau(signal+noise,gamma)
#' dy = diff(sdata)
#' ddy = diff(sdata,differences=2)
#' model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.1)
#' breaks = est.pair(vall=model2$vall,peak=model2$peak,gamma=gamma)$cp
#' slope = est.slope(x=(signal+noise),breaks=breaks)
cpTest = function(x,order,alpha,gamma,sigma,breaks,slope,untest,nu,is.constant,margin){
  if(!order %in% c(1,2))
    stop("Please specify 1 or 2 for the first or second derivative respectively")
  if(missing(gamma)) stop("gamma is required")
  if(missing(untest)) untest = NULL
  if(missing(is.constant)) is.constant = FALSE
  if(missing(nu)) nu = 0
  if(missing(margin)) margin = length(x)
  if(missing(sigma)) xi = sqrt(gamma^2+nu^2)
  else xi = sigma
  sigma = ifelse(order==1,sqrt(1/(4*xi^3*sqrt(pi))),sqrt(3/(8*xi^5*sqrt(pi))))
  kappa = ifelse(order==1,3/sqrt(5),15/sqrt(105))
  len = length(x)
  margins = unlist(lapply((0:(len %/% margin))*margin, function(t)(t-ceiling(1.5*gamma)):(t+ceiling(1.5*gamma))))
  untests = unique(unlist(lapply(untest, function(t)(t-ceiling(1.5*gamma)):(t+ceiling(1.5*gamma)))))
  # peak height distribution
  f1 = function(x) sqrt(3-kappa^2)/(sigma*sqrt(3))*dnorm(sqrt(3)*x/(sigma*sqrt(3-kappa^2)))+
    sqrt(2*pi)*kappa*x/(sigma^2*sqrt(3))*dnorm(x/sigma)*pnorm(kappa*x/(sigma*sqrt(3-kappa^2)))
  # minus slope
  f2 = function(locind){
    ty = vector()
    for(i in 1:length(slope)){
      if(i==1) seg = locind[locind<breaks[i]]
      else if(i==length(slope)) seg = locind[locind>=breaks[i-1]]
      else seg = locind[locind<breaks[i] & locind>=breaks[i-1]]
      ty = append(ty,x[seg]-slope[i])}
    return(ty)
  }
  #standardize: improve the accuracy of numerically computed p-value
  x = x/sigma
  # local extreme
  lmax = which.peaks(x)
  lmin = which.peaks(x,decreasing=T)
  if(order==1){
    if(is.constant) Ty = c(x[lmax],-x[lmin])
    else{
      if(missing(breaks) | missing(slope))
        stop("breaks and slope are required")
      if(length(breaks)+1 != length(slope))
        stop("breaks and slope should be matched")
      slope = slope/sigma  #standardize
      Ty = c(f2(lmax),-f2(lmin))
    }
  }
  else{
    lmax = setdiff(lmax,untests)
    lmin = setdiff(lmin,untests)
    Ty = c(x[lmax],-x[lmin])
  }
  # multiple testing
  sigma = 1 #standardize
  pval = unlist(lapply(Ty,function(x) integrate(f1,lower=x,upper=Inf)$value))
  pthresh = fdrBH(pval,alpha)
  if(pthresh==0){
    peak = vall = NULL; pthresh = 0
    warning("threshold for p-value is 0",call.=FALSE)
  }
  else{
    peak = intersect(c(lmax,lmin)[which(pval<=pthresh)],lmax)
    vall = intersect(c(lmax,lmin)[which(pval<=pthresh)],lmin)
    peak = peak[!(peak %in% margins)] # delete the margin points
    vall = vall[!(vall %in% margins)]
    if(length(peak)==0) peak = NULL
    if(length(vall)==0) vall = NULL
  }
  return(list(peak=peak,vall=vall,thresh=pthresh))
}
#' Detection of change points based on 'dSTEM' algorithm
#'
#' @param data vector of data sequence
#' @param type "I" if the change points are piecewise linear and continuous;
#'             "II-step" if the change points are piecewise constant and noncontinuous;
#'             "II-linear" if the change points are piecewise linear and noncontinuous;
#'             "mixture" if both type I and type II change points are include in \code{data}
#' @inheritParams cpTest
#' @inheritParams cpTest

#' @seealso \code{\link{cpTest}}
#' @return if type is 'mixture', the output is a list of type I and type II change points,
#'         otherwise, it is a list of change points

#' @export
#' @examples
#' ## piecewise linear signal
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = rep(0,7)
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5)/50
#' beta1 = c(beta1,-sum(beta1*(c(h[1],diff(h))))/(l-tail(h,1)))
#' signal = gen.signal(l,h,jump,beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' model = dstem(signal + noise,"I",gamma=gamma,alpha=0.05)
#' ## piecewise constant
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)
#' beta1 = rep(0,length(h)+1)
#' signal = gen.signal(l,h,jump,beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' model = dstem(signal + noise, "II-step",gamma,alpha=0.05)
#' ## piecewise linear with jump
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)*3
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5,-0.5)/50
#' signal = gen.signal(l=l,h=h,jump=jump,b1=beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' model = dstem(signal + noise, "II-linear",gamma,alpha=0.05)
dstem = function(data,type = c("I","II-step","II-linear","mixture"), gamma =20, alpha=0.05){
  type = match.arg(type)
  dy = diff(smth.gau(data,gamma))
  ddy = diff(dy)
  if (type == "I") {
    est = cpTest(x=ddy,order=2,gamma=gamma,alpha=alpha)
    out = list(vall = est$vall,peak = est$peak)
  }
  else if (type == "II-step") {
    est = cpTest(x=dy,order=1,alpha=alpha,gamma=gamma,is.constant=T)
    out = list(vall = est$vall,peak = est$peak)
  }
  else if (type == "II-linear") {
    model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=2*alpha)
    breaks = est.pair(model2$vall,model2$peak,gamma)$cp
    slope = est.slope(data,breaks)
    est = cpTest(x=dy,order=1,alpha=alpha,gamma=gamma,breaks=breaks,slope=slope)
    out = list(vall = est$vall, peak = est$peak)
  }
  # mixture
  else {
    model2 = cpTest(x=ddy,order=2,gamma=gamma,alpha=2*alpha)
    breaks = est.pair(model2$vall,model2$peak,gamma)$cp
    if(length(breaks)==0) breaks = floor(length(data)/2)
    slope = est.slope(data,breaks)
    model1 = cpTest(x=dy,order=1,alpha=alpha,gamma=gamma,breaks=breaks,slope=slope)
    jump_seg = unique(c(sapply(c(model1$peak,model1$vall),function(x) floor(x-2*gamma):ceiling(x+2*gamma))))
    est = cpTest(x=ddy,order=2,gamma=gamma,alpha=alpha,untest=jump_seg)
    out = list(type1 = list(vall=est$vall,peak=est$peak), type2 = list(vall=model1$vall,peak=model1$peak))
  }
  return(out)
}
#' Compute SNR of a certain change point location
#'
#' @param order order of derivative of data
#' @inheritParams smth.gau
#' @param is.jump logical value indicating if the location to be calculated is a jump point
#' @param jump jump height
#' @param diffb difference of the slopes on left and right sides of the location
#' @param addb sum of the slopes, only used when order is 1
#'
#' @return value of SNR
#' @export
snr = function(order, gamma, is.jump, jump, diffb, addb){
  if(!order %in% c(1,2))
    stop("Please specify 1 or 2 for the first or second derivative respectively")
  if(missing(addb)) addb = NULL
  if(missing(diffb)) diffb = NULL
  if(!is.jump) jump = 0
  if(order==1 && !is.jump)
    stop("SNR does not exist")
  sigma = ifelse(order==1,sqrt(1/(4*gamma^3*sqrt(pi))),sqrt(3/(8*gamma^5*sqrt(pi))))
  if(order==1){
    if(is.null(addb)) stop("addb is required")
    SNR = (addb/2 + jump/(gamma*sqrt(2*pi)))/sigma
  }
  else{
    if(is.null(diffb)) stop("diffb is required")
    if(!is.jump) SNR = 2*gamma^(3/2)/(sqrt(3)*pi^(1/4))*diffb
    else SNR = dnorm(1)/gamma*c(jump/gamma+diffb, -jump/gamma+diffb)/sigma
  }
  return(abs(SNR))
}

#' Compute TPR and FPR
#'
#' @param uh numerical vector of estimated change point locations
#' @param th numerical vector of true change point locations
#' @param b  location tolerance, usually specified as the bandwidth \code{gamma}
#'
#' @return a dataframe of \code{FDR} (FPR) and \code{Power} (TPR)
#' @export
Fdr = function(uh, th, b){
  confband = function(th,b) unique(unlist(lapply(th,function(x){seq(ceiling(x-b),floor(x+b))})))
  if (length(uh)==0) {
    FDR = 0
    Power = 0}
  else{
    n.tp = sum(uh %in% confband(th,b))
    FDR = 1 - n.tp/length(uh)
    Power = min(n.tp/length(th),1)
  }
  return(data.frame(matrix(c(FDR,Power),nrow=1,dimnames=list(NULL,c("FDR","Power")))))
}

#' Plot data sequence, the first and second-order derivatives, and their local extrema
#'
#' @param x numerical vector of signal or signal-plus-noise data
#' @param order order of derivative of data
#' @param icd.noise logical value indicating if \code{x} includes noise
#' @param H optional, vector of change-point locations
#'
#' @return a plot
#' @import MASS
#' @export
#' @examples
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)*3
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5,-0.5)/50
#' signal = gen.signal(l,h,jump,beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' sdata = smth.gau(signal+noise,gamma)
#' dy = diff(sdata)
#' ddy = diff(sdata,differences=2)
#' cp.plt(signal,0,FALSE)
#' points(signal+noise,col="grey")
#' cp.plt(dy,1,H=h)
#' cp.plt(ddy,2,H=h)
#'
cp.plt = function(x,order,icd.noise,H){
  if(!order %in% c(0,1,2)) stop("Please specify 0, 1 or 2 for original, first or second derivative")
  if(missing(icd.noise)) icd.noise = TRUE
  if(!is.logical(icd.noise)) stop(" 'icd.noise' is a logical value ")
  if(order==0){
    title = ifelse(icd.noise,"Signal Plus Noise","Signal")
    plot(x,xlab="t",ylab="",main=title,col="blue")
  }
  else{
    title = ifelse(order==1,"first-order derivative","second-order derivative")
    ylab = ifelse(order==1,"Dy","D2y")
    lmax = which.peaks(x)
    lmin = which.peaks(x, decreasing = T)
    if(!missing(H)){
      plot(x,type="l",lwd=1.5, xaxt="n",col="blue",xlab="t",ylab="",main=title)
      abline(v=H,lwd=1.5,lty=2)
      legend("topleft",legend="true location",lty=2,lwd=1.5,bty="n")
      axis(1,at=H)
    }
    else plot(x,type="l",lwd=1.5,col="blue",xlab="t",ylab=ylab,main=title)
    if(order==2) abline(h=0,lwd=1,col="grey")
    points(lmax, x[lmax], pch = 16, col = "green")
    points(lmin, x[lmin], pch = 16, col = "red")
  }
}

# f.int = function(x,b){
#   # Generate intervals for the detected change-point
#   # x is the estimate of change points
#   # b is the location torlerance
#   do.call(rbind,lapply(x,function(i) data.frame(starts=max(c(0,i-b)),ends=i+b)))
# }
# draw.rects <- function(interval, yrange, density = 10, col = "red") {
#   # Draw intervals of significance, as shaded rectangular areas, on the current plot.
#   # interval - a matrix or dataframe, interval of change points.
#   # yrange - vector of length two specifying the (lower, upper) vertical limit of the rectangles.
#   # density - density of the shading; try using 10 or 20.
#   # col - colour of the shading.
#   d = dim(interval)
#   for (i in 1:d[1]) {
#     rect(interval[i,1], yrange[1], interval[i,2],
#          yrange[2], density=density, col=col)
#   }
# }

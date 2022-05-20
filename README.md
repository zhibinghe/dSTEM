# R Package: dSTEM (differential Smoothing and TEsting of Maxima/Minima)

A new generic and fast approach for multiple change points detection in linear models. The method is designed for the data sequence 
modeled as piecewise linear signal plus noise (a stationary Gaussian process). 

# Toy Example of the Method

The core idea of our method is to transform the change points into local maxima/minima of the differential smoothed data. By performing 
multiple testing on the local maxima/minima, the significant ones are detected as change points.

For example, a continuous change point will be transformed into a local maximum/minimum in the second derivative of the smoothed data, 
while a noncontinuous change point will be transformed into a local maximum/minimum in the first derivative of the smoothed data.
A toy example is illustrated in [toy1.pdf](https://github.com/zhibinghe/ChangePoint/files/8721035/illu1.pdf) and [toy2.pdf](https://github.com/zhibinghe/ChangePoint/files/8721085/Illu2_1.pdf).

# Usage

The main function is 'cpTest()', which detect the change points in piecewise constant and piecewise linear signals.

``` r
#' Detection of change points based on 'dSTEM' algorithm
#'
#' @param data vector of data sequence
#' @param type "I" if the change points are piecewise linear and continuous;
#'             "II-step" if the change points are piecewise constant and noncontinuous;
#'             "II-linear" if the change points are piecewise linear and noncontinuous;
#'             "mixture" if both type I and type II change points are inclued in \code{data}
#' @inheritParams cpTest
#' @inheritParams cpTest

#' @seealso \code{\link{cpTest}}, \link[not]{features}
#' @return if type is 'mixture', the output is a list of type I and type II change points, 
#'         otherwise, it is a list of positive and negative change points

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
#'
#' ## piecewise constant
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)
#' beta1 = rep(0,length(h)+1)
#' signal = gen.signal(l,h,jump,beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' model = dstem(signal + noise, "II-step",gamma,alpha=0.05)
#'
#' ## piecewise linear with jump
#' l = 1200
#' h = seq(150,by=150,length.out=6)
#' jump = c(0,1.5,2,2.2,1.8,2,1.5)*3
#' beta1 = c(2,-1,2.5,-3,-0.2,2.5,-0.5)/50
#' signal = gen.signal(l=l,h=h,jump=jump,b1=beta1)
#' noise = rnorm(length(signal),0,1)
#' gamma = 25
#' model = dstem(signal + noise, "II-linear",gamma,alpha=0.05)
```

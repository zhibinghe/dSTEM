# R Package: dSTEM (differential Smoothing and TEsting of Maxima/Minima)

A new generic and fast approach for multiple change points detection in linear models. The method is designed for the data sequence 
modeled as piecewise linear signal plus noise (comes from a stationary Gaussian proccess). 

# Toy Example of the Method

The core idea of our method is to transform the change points into local maxima/minima of the differential smoothed data. By performing 
multiple testing on the local maxima/minima, the significant ones are detected as change points.

For example, a continuous change point will be transformed into a local maximum/minimum in the second derivative of the smoothed data, 
while a noncontinuous change point will be transformed into a local maximum/minimum in the first derivative of the smoothed data.
A toy example is illustrated in [toy1.pdf](https://github.com/zhibinghe/ChangePoint/files/8721035/illu1.pdf) and [toy2.pdf](https://github.com/zhibinghe/ChangePoint/files/8721085/Illu2_1.pdf).

# Usage

The main function is 'cpTest()', which detect the change points in piecewise constant and piecewise linear signals.

``` r
#' Mutilple testing for change-point based on 'dSTEM' algorithm
#'
#' @param x vector of data to be tested
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

#' @return a list of estimated change-point locations and threshold for p-value

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
model = cpTest(x=ddy,order=2,gamma=gamma,alpha=0.05)
```

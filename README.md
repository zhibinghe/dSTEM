# R Package: dSTEM (differential Smoothing and TEsting of Maxima/Minima)

A new generic and fast approach for multiple change points detection in linear models. The method is designed for the data sequence 
modeled as piecewise linear signal plus noise (comes from a stationary Gaussian proccess). 

# Toy Example of the Method

The core idea of our method is to transform the change points into local maxima/minima of the differential smoothed data. By performing 
multiple testing on the local maxima/minima, the significant ones are detected as change points.

For example, a continuous change point will be transformed into a local maximum/minimum in the second derivative of the smoothed data, 
while a noncontinuous change point will be transformed into a local maximum/minimum in the first derivative of the smoothed data.
A toy example is illustrated in [illu1.pdf](https://github.com/zhibinghe/ChangePoint/files/8721035/illu1.pdf).

# Usage


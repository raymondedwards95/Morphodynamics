function [Sk, As] = sk_as(eta)
% [Sk, As] = sk_as(eta)
% Computes the wave skewness Sk and wave asymmetry As from a detrended time-series of free-surface elevation eta
% Input:  eta   detrended time-series of free-surface elevation (m)
% Output: Sk    skewness
%         As    asymmetry
%

val_2 = power(mean(power(eta, 2)), 3/2);

Sk = mean(power(eta, 3)) / val_2;

As = mean(power(imag(hilbert(eta)), 3)) / val_2;

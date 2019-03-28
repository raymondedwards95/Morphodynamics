function [rms_value] = rms(x)
% Compute the root mean square of an array x
% Input:  x array
% Output: rootmeansquare

rms_value = sqrt(mean(power(x, 2)));

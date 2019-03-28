function [a] = acceleration(u, t)
% Computes the acceleration a from given timeseries of u (and t).
% Input:  u time-series of horizontal orbital velocity (m/s)
%         t corresponding time
% Output: a time-series of horizontal orbital acceleration (m/s2)

a = gradient(u, t);

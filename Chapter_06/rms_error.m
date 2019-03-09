function [err] = rms_error(X1, X2)
% Compute the root mean square error between two vectors X1 and X2
% returns \sqrt(1/N \sum^N_{i=1} ((X1_{i} - X2_{i})^2))
% input:  X1 - vector
%         X2 - vector
% output: err - root mean square error

err = sqrt(mean(power(X1 - X2, 2)));

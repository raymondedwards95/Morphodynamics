function [err] = rms_error(X1, X2)
% Compute the root mean square error between two vectors X1 and X2
% returns sqrt(1/N sum^N_{i=1} ((X1_{i} - X2_{i})^2))

%%% N = length(X2);
%%% square_diff = power(X1 - X2, 2);
%%% err = sqrt(1/N * sum(square_diff));

err = sqrt(mean(power(X1 - X2, 2)));

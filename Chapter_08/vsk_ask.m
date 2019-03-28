function [R, B] = vsk_ask(u, t)
% Computes the velocity skewness and the acceleration skewness parameters, R and beta in literature.
% Input:  u time-series of horizontal orbital velocity (m/s)
%         t corresponding time
% Output: R velocity skewness
%         B (beta) acceleration skewness

a = acceleration(u, t);

uc = abs(max(u));
ut = abs(min(u));
ac = abs(max(a));
at = abs(min(a));

R = uc ./ (uc + ut);
B = ac ./ (ac + at);

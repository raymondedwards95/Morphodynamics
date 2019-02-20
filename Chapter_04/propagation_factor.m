function [n] = propagation_factor(k, h)
% [n] = propagation_factor(k, h)
% Calculates the ratio n between the phase velocity and the group velocity
% input  k wave number (m-1)
%        h water depth (m)
% output n propagation factor ()

n = 1/2 * (1 + (2*k*h) / (sinh(2*k*h)));

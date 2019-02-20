function [cg] = group_velocity(T, h)
% Calculates the group velocity cg from the wave period T and the water depth h, 
% using cg = c * n, where c is the phase velocity and n is the propagation factor
% input  T wave period (s)
%        h water depth (m)
% output cg group velocity (m/s)

c = phase_velocity(T, h); % phase velocity
k = wave_number(T, h); % wave number, for n
n = propagation_factor(k, h); % propagation factor
cg = n * c;

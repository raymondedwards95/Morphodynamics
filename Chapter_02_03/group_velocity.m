function [cg] = group_velocity(T, h)
% Calculates the group velocity cg from the wave period T and the water depth h, 
% using cg = c * n, where c is the phase velocity and n is the propagation factor
% input  T wave period (s)
%        h water depth (m)
% output cg group velocity (m/s)

k = wave_number(T, h);
c = phase_velocity(T, h);
n = propagation_factor(k, h);
cg = n * c;

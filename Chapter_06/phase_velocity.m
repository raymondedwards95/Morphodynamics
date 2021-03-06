function [c] = phase_velocity(T, h)
% Calculates the phase velocity c from the wave period T and the water depth h, 
% from c = omega / k, where omega is the wave frequency
% input  T wave period (s)
%        h water depth (m)
% output c phase velocity (m/s)

omega = 2 * pi / T; % wave frequency
k = wave_number(T, h); % wave number
c = omega / k; % phase velocity

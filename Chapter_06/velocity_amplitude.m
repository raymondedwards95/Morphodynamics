function [Uw] = velocity_amplitude(h, Hrms, T)
% [Uw] = velocity_amplitude(h, Hrms, T)
% Computes the skewness and asymmetry of a wave using a parameterization defined in Reussink et al. (2012)
% Input:  h     water depth (m)
%         Hrms  root mean square wave height (m)
%         T     period (s)
% Output: Uw    velocity amplitude (m/s)
%

k = wave_number(T, h);

Uw = pi ./ T .* Hrms ./ sinh(k.*h);


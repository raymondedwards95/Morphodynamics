function [Ux] = undertow_magnitude(E, Er, c, h, rho, Hrms)
% Compute the magnitude of the undertow
% Input:  E     wave energy [J]
%         Er    roller energy [J]
%         c     phase velocity [m/s]
%         h     total water depth [m]
%         rho   water density [kg/m3]
%         Hrms  root mean square wave height [m]
% Output: Ux    magnitude of undertow [m/s]

M = (E + 2*Er) ./ c;
h_trough = h - Hrms ./ 2;

Ux = M ./ (rho * h_trough);

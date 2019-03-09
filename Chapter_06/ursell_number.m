function [Ur] = ursell_number(k, h, Hrms)
% [Ur] = ursell_number(k, h, Hrms)
% Computes the Ursell number Ur = 3ak/(4(kh)^3) from k, h and Hrms
% Input:  k     wave number (1/m)
%         h     water depth (m)
%         Hrms  root mean square wave height (m)
% Output: ur    Ursell number
%
%

a = 0.5 * sqrt(2) * Hrms;

Ur = 3/4 * a .* k ./ power(k.*h, 3);

function [k] = wave_number(T, h)
% [k] = wave_number(T, h)
% Calculates the wave number k from the wave period T and water depth h, using the approximation method of Gua (2002)
% input  T wave period (s)
%        h water depth (m)
% output k wave number (m-1)

g = 9.81; % gravitational acceleration
b = 2.4908; % approximation parameter

omega = 2 * pi / T; % wave frequency
x = h * omega / sqrt(g * h); % approximation 1
y = power(x, 2.) * power(1. - exp(-power(x, b)), -1./b); % approximation 2
k = y / h; % wave number

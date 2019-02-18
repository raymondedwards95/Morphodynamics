function [Hrms] = rms_height(H)
% [Hrms] = rms_height(H)
% Calculation of the root mean square waveheight H_rms of a list H with waveheights
% input  H    waveheights (m)
% output Hrms root mean square waveheight (m)

n = length(H); % number of elements
H_square = power(H, 2.); % take the square of all elements
Hrms = sqrt(1/n * sum(H_square)); % sum all squares, devide by n and take square-root

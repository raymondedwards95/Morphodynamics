function [Sk, As] = sk_as_ursell(Ur)
% [Sk, As] = sk_as_ursell(Ur)
% Computes the skewness and asymmetry of a wave using a parameterization defined in Reussink et al. (2012)
% Input:  Ur    Ursell number
% Output: Sk    skewness
%         As    asymmetry
%
a = (-0.471 - log10(Ur)) / 0.297;
B = 0.857 ./ (1 + exp(a));

Psi = -90 + 90 * tanh(0.815 ./ power(Ur, 0.672));
Psi = deg2rad(Psi);

As = B .* sin(Psi);
Sk = B .* cos(Psi);

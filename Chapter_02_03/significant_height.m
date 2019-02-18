function [Hs] = significant_height(H)
% [Hs] = significant_height(H)
% Calculation of the significant waveheight H_1/3 of a list H with waveheights
% input  H  waveheights (m)
% output Hs significant waveheight (m)

n = length(H); % number of elements
H_sorted = sort(H); % sort list
Hs = 0; % start with 0
for i=round(2*n/3)+1:1:n % for loop over the last 1/3 of the data
    Hs = Hs + H_sorted(i) / (1/3*n); % sum and divide by (1/3)*n
end

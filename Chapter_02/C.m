clear all
close all

f_s = 2;
f_N = f_s / 2;

fmin = [0, 0.005, 0.05];
fmax = [f_N, 0.05, f_N];

blocks = 15;

lowtide = load('lowTide.txt');
midtide = load('midTide.txt');
hightide = load('highTide.txt');

data = cat(3, lowtide, midtide, hightide);
clear lowtide midtide hightide

[n_s, n_p, n_t] = size(data);
assert(length(fmin) == length(fmax), 'fmin and fmax have different lengths');
n_f = min([length(fmin), length(fmax)]);

nfft = round(n_s ./ ((blocks + 1) / 2));

% Hm0_all = zeros(n_p, n_t, n_f);
m0_all = zeros(n_p, n_t, n_f);
H13 = zeros(n_p, n_t);
fp = zeros(n_p, n_t);

for i=1:1:n_p % loop over all positions
    for j=1:1:n_t % loop over all tides
        temp_data = data(:,i,j);
        [wave] = zero_crossing(temp_data, f_s);
        [S, f, e, c] = VarianceDensitySpectrum(temp_data, nfft, f_s);
        for k=1:1:n_f % loop over all pairs of fmin and fmax
            m0_all(i,j,k) = spectral_moment(f, S, fmin(k), fmax(k), 0);
        end
        H13(i,j) = significant_height(wave);

        [m, mi] = max(S);
        fp(i,j) = f(mi);
    end
end

Tp = 1 ./ fp;
Hm0_all = 4 * sqrt(m0_all);

Hm0 = Hm0_all(:,:,1);

figure;
scatter(reshape(Hm0, n_t*n_p, 1), reshape(H13, n_t*n_p, 1))
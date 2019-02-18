% Spectral analysis: Computation of spectral wave characteristics
% Chapter 2.2
% Script 1
%
% Computing spectral wave heights (characteristics) using zeroth-order moment
clear all
close all


%%% SETTINGS
f_s = 2; % sampling frequency
f_N = f_s / 2; % Nyquist frequency
fmin = [0, 0.005, 0.05]; % minimum frequencies for each range
fmax = [f_N, 0.05, f_N]; % maximum frequencies for each range
flabels = {'m_{0}'; 'm_{0,inf}'; 'm_{0,ss}'}; % names for each range

blocks = 15; % number of blocks


%%% READ DATA
lowtide = load('lowTide.txt');
midtide = load('midTide.txt');
hightide = load('highTide.txt');

data = cat(3, lowtide, midtide, hightide); % stack data in a new dimension, result contains 3 axes
clear lowtide midtide hightide % remove unused variables

bedprofile = load('prof1018.txt');
x_b = bedprofile(:,1); % x-coordinates bed
y_b = bedprofile(:,2); % y-coordinates bed
x_p = [4478, 4765, 4790, 4814, 4835]; % x-coordinates positions
clear bedprofile


%%% PREPARE CALCULATIONS
[n_s, n_p, n_t] = size(data); % number of samples, number of positions, number of datasets (tides)

assert(length(fmin) == length(fmax), 'fmin and fmax have different lengths'); % check if fmin and fmax have the same length
n_f = min([length(fmin), length(fmax)]); % if not, use the minimum length

nfft = round(n_s ./ ((blocks + 1) / 2)); % number of elements in each block

m0_all = zeros(n_p, n_t, n_f); % save for each position, for each tide and for each frequency range
H13 = zeros(n_p, n_t); % save for all positions and for all tides
fp = zeros(n_p, n_t); % all positions, all tides


%%% CALCULATIONS
for i = 1:n_p % loop over all positions
    for j = 1:n_t % loop over all tides
        temp_data = data(:,i,j); % take subset of data for calculations

        [S, f, ~, ~] = VarianceDensitySpectrum(temp_data, nfft, f_s);
        for k = 1:n_f % loop over all pairs of fmin and fmax
            m0_all(i,j,k) = spectral_moment(f, S, fmin(k), fmax(k), 0);
        end

        [wave] = zero_crossing(temp_data, f_s); % find the waves to calculate the H13
        H13(i,j) = significant_height(wave);

        [~, mi] = max(S); % find maximum and the corresponding index
        fp(i,j) = f(mi); % find the corresponding frequency
    end
end
clear mi S f temp_data i j k wave % remove unused variables

Tp = 2 * pi ./ fp; % T = 2 pi / f % elementwise matrix operation
Hm0_all = 4 * sqrt(m0_all); % Hm0 = 4*sqrt(m0)


%%% EXTRACT RESULTS
Hm0 = Hm0_all(:,:,1);
Hm0_inf = Hm0_all(:,:,2);
Hm0_ss = Hm0_all(:,:,3);


%%% FIGURE 1
% comparison of Hm0 and H13
Hm0f = reshape(Hm0, n_t*n_p, 1);
H13f = reshape(H13, n_t*n_p, 1);
coeff = polyfit(Hm0f, H13f, 1);
figure;
hold on
scatter(Hm0f, H13f)
plot(Hm0f, polyval(coeff, Hm0f))
hold off
xlabel('H_{m0}')
ylabel('H_{1/3}')
legend('data', ['fit: slope = ', num2str(coeff(1))], 'Location', 'NorthWest')
title('H_{1/3} and H_{m0}')
saveas(gcf, 'figures/2_2_spectral_compare.png')


%%% FIGURE 2
% cross-shore evolution for all frequency ranges
tide = 1; % tide number: 1=low, 2=mid, 3=high
ylims = [2, 0.5, 2]; % y ranges from 0 to *. Length of array should be the same as length of fmin, fmax, flabels
x_left = 4350; % minimum value for x in plots
x_right = max(x_p)+10; % maximum value

figure;
sgtitle('Spectral wave heights') % title of figure

for k = 1:n_f
    subplot(n_f+1, 1, k);
    hold on
    scatter(x_p, Hm0_all(:,tide,k), '+') % plot data as points
    % plot(x_p, Hm0_all(:,tide,k), '--') % connect points with straight lines
    hold off
    xlim([x_left, x_right])
    ylim([0, ylims(k)])
    set(gca, 'xticklabel', []) % remove xticklabels
    ylabel([flabels(k), ' [m]'])
    legend([flabels(k)], 'Location', 'NorthWest')
end

subplot(n_f+1, 1, n_f+1)
hold on
plot(x_b, y_b) % bed profile
scatter(x_p, interp1(x_b, y_b, x_p)) % dots, using interp to get y-value for x of all positions
hold off
xlim([x_left, x_right])
ylim([-7, 0])
title('Bed profile')
xlabel('Distance from bouy [m]')
ylabel('Elevation [m]')
saveas(gcf, 'figures/2_2_spectral_low.png')

save('data_spectral', 'Hm0', 'Hm0_ss', 'Hm0_inf', 'Tp', 'H13')
clear k ylims x_left x_right % remove not-important variables


%%% FIGURE 3

% Linear wave theory: Egmond dataset
% Chapter 3.2
% Script 1
%
% Wavelenghts derived from the Egmond dataset
clear all
close all


%% DATA
T_13 = [7.58, 6.69, 5.54]; % periods as given in table 9.1
h = load('MeanWaterDepth.txt'); % load mean water depth during low, mid and high tides at different positions
[n_p, n_t] = size(h); % number of positions, number of tides

bedprofile = load('prof1018.txt'); % load bed profile
x_b = bedprofile(:,1); % x-coordinates bed
y_b = bedprofile(:,2); % y-coordinates bed
x_p = [4478, 4765, 4790, 4814, 4835]; % x-coordinates positions
clear bedprofile


%% PREPARE CALCULATIONS
L = zeros(n_p, n_t); % pre-allocate array


%% CALCULATIONS
for i = 1:n_t % loop over all tides
    for j = 1:n_p % loop over all positions
        L(j,i) = 2 * pi / wave_number(T_13(i), h(j,i)); % calculate wavelength
    end
end

hL = h ./ L; % calculate the ratio h/L (elementwise)


%% FIGURE
x_left = 4450; % minimum value for x in plots
x_right = 4950; % max(x_p)+100; % maximum value

figure;
pbaspect([2 1 1]) % set aspect ratio (?)

subplot(3,1,1);
hold on
for i = 1:n_t % loop over all tides
    scatter(x_p, L(:,i), '+') % plot L (points)
end
hold off
xlim([x_left, x_right])
ylim([20, 45])
set(gca, 'xticklabel', []) % remove xticklabels
title('Wavelength')
ylabel('L [m]')
legend('low', 'mid', 'high')

subplot(3,1,2);
hold on
for i = 1:n_t % loop over all tides
    scatter(x_p, hL(:,i), '+') % plot h/L (points)
end
plot([x_left, x_right], 0.05*ones(2, 1), '--k') % plot horizontal line for shallow water
plot([x_left, x_right], 0.5*ones(2, 1), '--k') % plot horizontal line for deep water
hold off
xlim([x_left, x_right])
ylim([0, 0.2])
title('Ratio h/L')
ylabel('h/L')
set(gca, 'xticklabel', []) % remove xticklabels
legend('low', 'mid', 'high')

subplot(3,1,3)
hold on
plot(x_b, y_b) % bed profile
scatter(x_p, interp1(x_b, y_b, x_p)) % dots, using interp to get y-value for x of all positions
hold off
xlim([x_left, x_right])
ylim([-7, 0])
title('Bed profile')
xlabel('Distance from bouy [m]')
ylabel('z [m]')

saveas(gcf, 'figures/3_2_wavelengths.png')

clear i j x_b x_left x_p x_right y_b

% Linear wave theory: Egmond dataset
% Chapter 3.2
% Script 1
%
% 
clear all
close all


%%% SETTINGS


%%% DATA
T_13 = [7.58, 6.69, 5.54];
h = load('MeanWaterDepth.txt');
[n_p, n_t] = size(h);

bedprofile = load('prof1018.txt');
x_b = bedprofile(:,1); % x-coordinates bed
y_b = bedprofile(:,2); % y-coordinates bed
x_p = [4478, 4765, 4790, 4814, 4835]; % x-coordinates positions
clear bedprofile


%%% PREPARE CALCULATIONS
L = zeros(n_p, n_t);


%%% CALCULATIONS
for i = 1:n_t
    for j = 1:n_p
        L(j,i) = 2 * pi / wave_number(T_13(i), h(j,i));
    end
end

hL = h ./ L;


%%% FIGURE
x_left = 4450; % minimum value for x in plots
x_right = max(x_p)+10; % maximum value

figure;

subplot(3,1,1);
hold on
for i = 1:n_t
    scatter(x_p, L(:,i), '+') % plot data as points
end
hold off
xlim([x_left, x_right])
set(gca, 'xticklabel', []) % remove xticklabels
ylabel('L [m]')

subplot(3,1,2);
hold on
for i = 1:n_t
    scatter(x_p, hL(:,i), '+') % plot data as points
end
plot([x_left, x_right], 0.05*ones(2, 1), 'k')
plot([x_left, x_right], 0.5*ones(2, 1), 'k')
hold off
xlim([x_left, x_right])
ylim([0, 0.2])
ylabel('h/L')
set(gca, 'xticklabel', []) % remove xticklabels

subplot(3,1,3)
hold on
plot(x_b, y_b) % bed profile
scatter(x_p, interp1(x_b, y_b, x_p)) % dots, using interp to get y-value for x of all positions
hold off
xlim([x_left, x_right])
ylim([-7, 0])
title('Bed profile')
xlabel('Distance from bouy [m]')
ylabel('H [m]')
% Time-series analysis: pre-processing
% Chapter 1.1
% Script 1
%
% Processing of sensor data
% Plotting of data
close all
clear all
disp('### Chapter 1.1: script 1 ###')


%%% SETTINGS
sampling = 4 % Hz = 1/s % sampling rate of sensor
timestep = 1/sampling % s % 1 divided by sampling rate
a1 = 19.97370 % Pa/mV % calibration parameter % p = ax + b
b1 = -49.95959 % Pa % calibration parameter
hb1 = 1.45 % m % height of sensor


%%% CONSTANTS
rho = 1025 % kg/m3 % water density
g = 9.81 % m/s2 # gravitational acceleration


%%% READ DATA
calibp1 = load('calibP1.txt'); % read data
number = length(calibp1) % number of elements in data
duration = number / sampling / 60 % duration in minutes


%%% CONVERT DATA
p1 = a1 * calibp1 + b1; % pressure signal
ha1 = p1 / rho / g; % height of water column above
h1 = ha1 + hb1; % total water depth


%%% PROCESS DATA 
h1_mean = mean(h1) % average water depth at the sensor
h1_d = detrend(h1); % h1 without trend

t1 = (timestep:timestep:duration*60) - timestep; % time starting at 0 and specified timestep


%%% PLOT DATA
disp('Creating figures')
% plot water depth
figure;
title('Water depth')
xlabel('Time [s]');
ylabel('Depth [m]');
xlim([min(t1), max(t1)])
ylim([2, 6])

hold on
plot(t1, h1);
plot([min(t1), max(t1)], [h1_mean, h1_mean], 'LineWidth', 2) % plot horizontal line
line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k') % plot y=0
hold off
legend('P1', 'P1 - average water depth')
saveas(gcf, 'figures/1_1_waterdepth.png')

% plot water depth variations
figure;
subplot(211)
plot(t1, h1-h1_mean) % plot deviations from the mean
line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k') % plot y=0
legend('P1')
xlim([min(t1), max(t1)])
ylim([-2, 2])
xlabel('Time [s]')
ylabel('Depth deviation [m]')
title('Depth variations')

subplot(212)
plot(t1, h1_d) % plot deviations from a linear fit
line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k') % plot y=0
legend('P1')
xlim([min(t1), max(t1)])
ylim([-2, 2])
xlabel('Time [s]')
ylabel('Depth deviation [m]')
title('Detrended depth')
saveas(gcf, 'figures/1_1_detrend.png')

disp('### End ###')

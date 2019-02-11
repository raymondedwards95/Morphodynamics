% Time-series analysis: wave statistics - preliminary visualisations
% Chapter 1.2
% Script 1
%
% Loading multiple sources of data
% Plotting subsets of data
close all
clear all
disp('### Chapter 1.2: script 1 ###')


%%% SETTINGS
sampling = 2 % Hz
timestep = 1/sampling


%%% READ DATA
low = load('lowTide.txt');
high = load('highTide.txt');

[ml nl] = size(low);
[mh nh] = size(high);

tl = (0:1:ml-1) * timestep;
th = (0:1:mh-1) * timestep;


%%% PLOT DATA
disp('Creating figures')
figure;
subplot(321)
plot(tl, low(:,1))
title('P1 - low tide')
xlabel('Time [s]')
ylabel('SSE [m]')
xlim([min(tl), max(tl)])
ylim([-2, 2])

subplot(323)
plot(tl, low(:,2))
title('P3 - low tide')
xlabel('Time [s]')
ylabel('SSE [m]')
xlim([min(tl), max(tl)])
ylim([-1.5, 1.5])

subplot(325)
plot(tl, low(:,5))
title('P6 - low tide')
xlabel('Time [s]')
ylabel('SSE [m]')
xlim([min(tl), max(tl)])
ylim([-1, 1])

subplot(322)
plot(th, high(:,1))
title('P1 - high tide')
xlabel('Time [s]')
ylabel('SSE [m]')
xlim([min(th), max(th)])
ylim([-2, 2])

subplot(324)
plot(th, high(:,2))
title('P3 - high tide')
xlabel('Time [s]')
ylabel('SSE [m]')
xlim([min(th), max(th)])
ylim([-1.5, 1.5])

subplot(326)
plot(th, high(:,5))
title('P6 - high tide')
xlabel('Time [s]')
ylabel('SSE [m]')
xlim([min(th), max(th)])
ylim([-1, 1])

saveas(gcf, 'figures/1_2_sse_multi.png')

disp('### End ###')

% Non-linear wave transformation: Observations of skewness and asymmetry
% Chapter 6.1
% Script 1
%
close all
clear all


%% SETTINGS
n_t = 3; % number of tides
n_p = 5; % number of locations

% sensors
xp16 = [4478, 4765, 4790, 4814, 4835]; % location of sensors 1,3-6

% variables
Ur_th = logspace(-2, 2, 100); % create logaritmic vector from 10^(-2) to 10^(2)


%% LOAD DATA
eta = cat(3, ...
          load('lowTide.txt'), ...
          load('midTide.txt'), ...
          load('highTide.txt') ...
         ); % free surface elevation
h = load('MeanWaterDepth.txt'); % water depth
param_egmond = load('data/52_ParametersEgmond.mat'); % values from table 9.1
T13 = param_egmond.T13; % significant period
Hdata = load('data/12_StatisticsEgmond'); % observed wave height (chapter 1.2)
Hrms = Hdata.Hrms_tot; % root mean square wave height

clear param_egmond Hdata


%% PREPARATIONS
% pre-allocation
As_eg = zeros(n_p, n_t);
Sk_eg = zeros(n_p, n_t);
Ur_eg = zeros(n_p, n_t);
k_eg = zeros(n_p, n_t);


%% CALCULATIONS
% theory
[Sk_th, As_th] = sk_as_ursell(Ur_th); % computes skewness and asymmetry from ursell number

% egmond
for i = 1:n_t % loop over all tides
    for j = 1:n_p % loop over all positions
        % compute wave number
        k_eg(j,i) = wave_number(T13(i), h(j,i));

        % compute skewness and asymmetry
        [Sk_eg(j,i), As_eg(j,i)] = sk_as(eta(:,j,i));

        % compute ursell number
        Ur_eg(j,i) = ursell_number(k_eg(j,i), h(j,i), Hrms(j,i));
    end
end


%% FIGURES
figure

% skewness
subplot(2,1,1)
semilogx(Ur_th, Sk_th)
hold on
grid on
box on
for i = 1:n_t
    scatter(Ur_eg(:,i), Sk_eg(:,i), '+')
end
ylim([-0.1, 1.2])
% xlabel('Ur')
ylabel('Sk')
title('Skewness')
legend('Fit', 'Low tide', 'Mid tide', 'High tide')
set(gca, 'xticklabel', []) % remove xticklabels

% asymmetry
subplot(2,1,2)
semilogx(Ur_th, As_th)
hold on
grid on
box on
for i = 1:n_t
    scatter(Ur_eg(:,i), As_eg(:,i), 'x')
end
ylim([-1, 0.1])
xlabel('Ur')
ylabel('As')
title('Asymmetry')
legend('Fit', 'Low tide', 'Mid tide', 'High tide')

% save fig
print('figures/61_theory_fit', '-dpng', '-r300')
print('figures/61_theory_fit', '-depsc', '-r300')


%% SAVE DATA
save('data/61_WaveShapesEg', 'Sk_eg', 'As_eg', 'Ur_eg')


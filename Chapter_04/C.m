% Modelling cross-shore wave transformation: Model/data comparison
% Chapter 4.3
% Script 1
%
close all
clear all


%% SETTINGS
% number of situations
n_t = 3; % low, mid, high

% off-shore conditions
Zeta = [-0.45, 0.09, 0.91]; % water levels during low tide, mid tide and high tide (m NAP)
T13 = [7.58, 6.69, 5.54]; % significant period (s)
H13 = [1.70, 2.25, 1.69]; % significant wave height (m)
theta = [-36, 39, 36]; % incident angle (deg)

% model parameter
hmin = 0.2; % Minimal water depth for computation (we stop the computation when h<hmin)

% locations of sensors
xp = [4478, 4765, 4790, 4814, 4835];

% visualisations
tide_names = {'Low tide', 'Mid tide', 'High tide'};
xlims = [4400, 5000];
ylims = [0, 1.6];


%% LOAD DATA
bedprofile = load('prof1018.txt');
xb = bedprofile(:,1);
zb = bedprofile(:,2);

H_obs = load('StatisticsEgmond.mat'); % observed data, Hrms values are in H_obs.Hrms_tot
Hrms_obs = H_obs.Hrms_tot;
clear H_obs


%% PREPARE CALCULATIONS
% Conversion of parameters
Hrms = H13 / sqrt(2); % H_{1/3} ~ sqrt{2} H_{rms}
T0 = T13; % T ~ T_{1/3}

% visualisations
if xlims(2) > max(xb)
    xlims = max(xb);
end

% pre-allocation
Hrms_model = zeros(5,3);
rms_err = zeros(1,3);


%% CALCULATIONS
for i = n_t:-1:1
    data(i) = BJmodel(Hrms(i), T0(i), Zeta(i), theta(i), bedprofile, hmin);

    Hrms_model(:,i) = interp1(data(i).x, data(i).Hrms, xp);

    rms_err(i) = rms_error(Hrms_obs(:,i), Hrms_model(:,i));
end

rms_err_tot = rms_error(Hrms_obs(:), Hrms_model(:));


%% VISUALISATIONS
figure
for i = 1:n_t
    subplot(4,1,i)
    hold on
    plot(data(i).x, data(i).Hrms) % model data
    scatter(xp, Hrms_obs(:,i), '+r') % observations
    hold off
    title(tide_names(i))
    xlim(xlims)
    ylim(ylims)
    ylabel('H_{rms} [m]')
    set(gca, 'xticklabel', []) % remove ticks
end


subplot(4,1,4)
hold on
for i = 1:n_t
    plot([min(xb), max(xb)], Zeta(i)*ones(2,1), '-.') % sea level
end
plot(xb, zb, 'k') % bed profile
scatter(xp, interp1(xb, zb, xp), 'r') % dots, using interp to get y-value for x of all positions
hold off
title('Bed profile')
legend(tide_names, 'Location', 'SouthEast')
xlim(xlims)
ylim([-7, 2])
xlabel('Cross-shore distance [m]')
ylabel('z [m]')

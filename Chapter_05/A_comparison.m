% Modelling of the alongshore flow: Data/model comparison for Egmond
% Chapter 5.2
% Script 1
%
clear all
close all


%% SETTINGS and PARAMETERS
n_t = 3;
tides = {'Low tide', 'Mid tide', 'High tide'};

Zeta = [-0.45, 0.09, 0.91]; % water levels during low tide, mid tide and high tide (m NAP)
H13 = [1.70, 2.25, 1.69]; % significant wave height (m)
T13 = [7.58, 6.69, 5.54]; % significant period (s)
theta = [-36, 39, 36]; % angle of wave incidence (degrees)
w_long = [-2.26, -3.27, 4.68]; % alongshore wind speed (m/s)
w_cross = [8.86, 8.55, 1.16]; % cross-shore wind speed (m/s)
dzeta_dy = [-7.33e-6, -1.76e-5, 8.57e-7]; % alongshore mean surface slope

% sensors
xp16 = [4478, 4765, 4790, 4814, 4835]; % location of sensors 1,3-6
xp38 = [4765, 4790, 4814, 4835, 4860, 4889]; % location of sensors 3-8

% constants
ka = 0.022; % apparent bed roughness (m)
nu = 0.5; % viscosity (m^2/s)

% model paramters
hmin = 0.2; % Minimal water depth for computation (we stop the computation when h<hmin)

% visualizations
xlims = [4200, 5200];
sz = 50; % scatter point size (default 36)
lw = 0.75; % plot line width (default 0.5)


%% LOAD DATA
bedprofile = load('prof1018.txt'); % load profile
v_obs = load('vEgmond.txt'); % load mean alongshore current data (at P3 to P8)
egmond_data = load('data/12_StatisticsEgmond'); % load observed wave height data
Hrms_obs = egmond_data.Hrms_tot;
clear egmond_data


%% PREPARATIONS
% convert parameters
Hrms = H13 / sqrt(2); % H_{1/3} ~ sqrt{2} H_{rms}
T0 = T13; % T ~ T_{1/3}
xb = bedprofile(:,1);
zb = bedprofile(:,2);

% pre-allocation
v_model = zeros(length(bedprofile), n_t);
v_model_interp = zeros(size(v_obs));

err_v = zeros(n_t+1, 1);

% visualizations
if xlims(2) > max(xb) % set xlim, depending on given value and bed profile
    xlims(2) = max(xb);
end


%% CALCULATIONS
% for i = n_t:-1:1 % loop over all tides, backwards for pre-allocation
for i = 1:n_t
    wavedata(i) = BJmodel(Hrms(i), T0(i), Zeta(i), theta(i), bedprofile, hmin); % calculate wave evolution

    v_model(:,i) = longshoreCurrent(bedprofile, dzeta_dy(i), w_cross(i), w_long(i), wavedata(i).c, wavedata(i).theta, wavedata(i).Dr, wavedata(i).ht, wavedata(i).st, ka, nu); % calculate alongshore current

    v_model_interp(:,i) = interp1(wavedata(i).x, v_model(:,i), xp38);

    err_v(i) = rms_error(v_obs(:,i), v_model_interp(:,i));
end
err_v(n_t+1) = rms_error(v_obs(:), v_model_interp(:));


%% VISUALIZATIONS
for i = 1:n_t
    figure
    sgtitle(tides(i))

    
    % Hrms
    subplot(3,1,1)
    box on
    grid on
    hold on
    plot(wavedata(i).x, wavedata(i).Hrms, 'LineWidth', lw) % show wave model
    scatter(xp16, Hrms_obs(:,i), sz, '+r') % show wave data points

    title('Root mean square wave height')
    legend('Model', 'Measurements', 'Location', 'SouthWest')
    xlim(xlims)
    ylim([0, 1.1*max(wavedata(i).Hrms)])
    ylabel('H_{rms} [m]')
    set(gca, 'xticklabel', []) % remove xticklabels


    % alongshore current
    subplot(3,1,2)
    box on
    grid on
    hold on
    plot(wavedata(i).x, v_model(:,i), 'LineWidth', lw) % show velocity model
    scatter(xp38, v_obs(:,i), sz, 'xr') % show velocity data points
    plot(xlims, zeros(2,1), '--k') % show v=0

    title('Alongshore current velocity')
    % legend('Model', 'Measurements', 'Location', 'NorthWest')
    xlim(xlims)
    % ylim([min(v_model(:,i)), max(v_model(:,i))])
    ylabel('v [m/s]')
    set(gca, 'xticklabel', []) % remove xticklabels


    % bed profile
    subplot(3,1,3)
    box on
    grid on
    hold on
    plot(xlims, Zeta(i)*ones(2,1), '-.') % show water level
    scatter(xp16, interp1(xb, zb, xp16), sz, '+r') % show wave data points
    scatter(xp38, interp1(xb, zb, xp38), sz, 'xr') % show velocity data points
    plot(xb, zb, 'k', 'LineWidth', lw)

    title('Bed profile')
    legend('Sea level', 'Wave measurement sensors', 'Current measurement sensors', 'Location', 'SouthEast')
    xlim(xlims)
    ylim([-8, 2])
    xlabel('x [m]')
    ylabel('z [m]')


end


%% SAVE DATA
save('data/52_ParametersEgmond', 'Zeta', 'H13', 'T13', 'theta', 'w_long', 'w_cross', 'dzeta_dy');
save('data/52_BJoutputAllTides', 'wavedata');
save('data/52_LongCurrentAllTides', 'v_model', 'xb', 'zb')


clear sz lw

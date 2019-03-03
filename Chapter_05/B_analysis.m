% Modelling of the alongshore flow: Complementary analysis
% Chapter 5.3
% Script 1
%
clear all
close all


%% SETTINGS and PARAMETERS
% number of situations:
n_t = 4;
% 1 normal
% 2 no wind
% 3 no wind, no tide
% 4 normal, angle=0

Zeta = [-0.45]; % water levels during low tide, mid tide and high tide (m NAP)
H13 = [1.70]; % significant wave height (m)
T13 = [7.58]; % significant period (s)
theta = [-36]; % angle of wave incidence (degrees)
w_long = [-2.26]; % alongshore wind speed (m/s)
w_cross = [8.86]; % cross-shore wind speed (m/s)
dzeta_dy = [-7.33e-6]; % alongshore mean surface slope

% repeat all matrices n_t times
dzeta_dy = repmat(dzeta_dy, 1, n_t);
w_cross = repmat(w_cross, 1, n_t);
w_long = repmat(w_long, 1, n_t);
theta = repmat(theta, 1, n_t);
Zeta = repmat(Zeta, 1, n_t);
T13 = repmat(T13, 1, n_t);
H13 = repmat(H13, 1, n_t);

% replace elements for different situations
w_cross(2:3) = [0,0]; % situations 2,3 without wind
w_long(2:3) = [0,0]; % situations 2,3 without wind
dzeta_dy(3) = [0]; % situation 3 without tides
theta(4) = [0]; % situation 4 with theta is zero

% sensors
xp16 = [4478, 4765, 4790, 4814, 4835]; % location of sensors 1,3-6
xp38 = [4765, 4790, 4814, 4835, 4860, 4889]; % location of sensors 3-8

% constants
ka = 0.022; % apparent bed roughness (m)
nu = 0.5; % viscosity (m2/s)

% model paramters
hmin = 0.2; % Minimal water depth for computation (we stop the computation when h<hmin)

% visualizations
xlims = [4000, 5200];
sz = 50; % scatter point size (default 36)
lw = 0.75; % plot line width (default 0.5)


%% READ DATA
bedprofile = load('prof1018.txt'); % load profile
v_obs = load('vEgmond.txt'); % load mean alongshore current data (at P3 to P8)


%% PREPARATIONS
% convert parameters
Hrms = H13 / sqrt(2); % H_{1/3} ~ sqrt{2} H_{rms}
T0 = T13; % T ~ T_{1/3}
xb = bedprofile(:,1);
zb = bedprofile(:,2);

% pre-allocation
v_model = zeros(length(bedprofile), n_t);

% visualizations
if xlims(2) > max(xb) % set xlim, depending on given value and bed profile
    xlims(2) = max(xb);
end


%% CALCULATIONS
for i = n_t:-1:1 % loop over all tides, backwards for pre-allocation
    wavedata(i) = BJmodel(Hrms(i), T0(i), Zeta(i), theta(i), bedprofile, hmin); % calculate wave evolution

    v_model(:,i) = longshoreCurrent(bedprofile, dzeta_dy(i), w_cross(i), w_long(i), wavedata(i).c, wavedata(i).theta, wavedata(i).Dr, wavedata(i).ht, wavedata(i).st, ka, nu); % calculate alongshore current
end


%% VISUALIZATIONS
figure
sgtitle('Low tide')

% current
subplot(2,1,1)
box on
grid on
hold on
for i=1:n_t
    plot(wavedata(i).x, v_model(:,i), 'LineWidth', lw)
end
scatter(xp38, v_obs(:,1), sz, 'xr') % show velocity data points
plot(xlims, zeros(2,1), '--k') % show v=0

title('Alongshore current velocity')
legend('All forcings', 'No wind', 'No wind, no tides', 'Normal incidence', 'Measurements', 'Location', 'SouthWest')
xlim(xlims)
ylabel('v [m/s]')
set(gca, 'xticklabel', []) % remove xticklabels


% bed profile
subplot(2,1,2)
box on
grid on
hold on
plot(xlims, Zeta(i)*ones(2,1), '-.') % show water level
% scatter(xp16, interp1(xb, zb, xp16), sz, '+r') % show wave data points
scatter(xp38, interp1(xb, zb, xp38), sz, 'xr') % show velocity data points
plot(xb, zb, 'k', 'LineWidth', lw)

title('Bed profile')
legend('Sea level', 'Sensors', 'Location', 'NorthWest')
xlim(xlims)
ylim([-8, 2])
xlabel('x [m]')
ylabel('z [m]')

% save figure
print('figures/5_3_compare', '-dpng')


%% SAVE DATA
% save('data/53_BJoutputLowTideMod', 'wavedata');
% save('data/53_LongCurrentLowTideMod', 'v_model', 'xb', 'zb')


clear sz lw

% Non-linear wave transformation: Modelling Sk and As
% Chapter 6.2
% Script 1
%
close all
clear all


%% SETTINGS
% tides
n_t = 3;
tides = {'Low tide', 'Mid tide', 'High tide'};

% model
hmin = 0.2; % Minimal water depth for computation (we stop the computation when h<hmin)

% sensors
xp16 = [4478, 4765, 4790, 4814, 4835]; % location of sensors 1,3-6

% visualizations
xlims = [4000, 5200];

%% LOAD DATA
load('data/52_ParametersEgmond.mat', 'H13', 'T13', 'Zeta', 'theta'); % load H13, T13, Zeta, theta
load('data/61_WaveShapesEg.mat', 'Sk_eg', 'As_eg', 'Ur_eg'); % load Sk_eg, As_eg, Ur_eg

bedprofile = load('prof1018.txt');
xb = bedprofile(:,1);
zb = bedprofile(:,2);


%% PREPARE
n_x = length(xb); % number of points
n_p16 = length(xp16); % number of sensors

% convert parameters
Hrms = H13 / sqrt(2); % H_{1/3} ~ sqrt{2} H_{rms}
T0 = T13; % T ~ T_{1/3}

% pre-allocation
Ur = zeros(n_x, n_t);
Sk = zeros(n_x, n_t);
As = zeros(n_x, n_t);

Sk_model = zeros(n_p16, n_t);
As_model = zeros(n_p16, n_t);

Sk_err = zeros(n_t+1, 1);
As_err = zeros(n_t+1, 1);

% visualizations
if xlims(2) > max(xb) % set xlim, depending on given value and bed profile
    xlims(2) = max(xb);
end


%% CALCULATIONS
for i = n_t:-1:1
    % cross-shore evolution
    wavedata(i) = BJmodel(Hrms(i), T0(i), Zeta(i), theta(i), bedprofile, hmin); % calculate wave evolution

    Ur(:,i) = ursell_number(wavedata(i).k, wavedata(i).ht, wavedata(i).Hrms); % calculate the ursell number for all points

    [Sk(:,i), As(:,i)] = sk_as_ursell(Ur(:,i)); % compute skewness and asymmetry from ursell

    % errors
    Sk_model(:,i) = interp1(xb, Sk(:,i), xp16);
    As_model(:,i) = interp1(xb, As(:,i), xp16);

    Sk_err(i) = rms_error(Sk_model(:,i), Sk_eg(:,i)); % compute error (rms)
    As_err(i) = rms_error(As_model(:,i), As_eg(:,i)); % compute error (rms)
end
Sk_err(n_t+1) = rms_error(Sk_model(:), Sk_eg(:));
As_err(n_t+1) = rms_error(As_model(:), As_eg(:));


%% FIGURE
for i = 1:n_t
    figure
    sgtitle(tides(i))

    % skewness
    subplot(3,1,1)
    box on
    grid on
    hold on
    plot(wavedata(i).x, Sk(:,i))
    scatter(xp16, Sk_eg(:,i), '+')

    title('Skewness')
    xlim(xlims)
    ylim([-0.1, 1.1])
    ylabel('Sk')
    set(gca, 'xticklabel', []) % remove xticklabels

    % asymmetry
    subplot(3,1,2)
    box on
    grid on
    hold on
    plot(wavedata(i).x, As(:,i))
    scatter(xp16, As_eg(:,i), 'x')

    title('Asymmetry')
    xlim(xlims)
    ylim([-1.1, 0.1])
    ylabel('As')
    set(gca, 'xticklabel', []) % remove xticklabels

    % bed profile
    subplot(3,1,3)
    box on
    grid on
    hold on
    plot(xlims, Zeta(i)*ones(2,1), '-.') % show water level
    scatter(xp16, interp1(xb, zb, xp16), '+r') % show wave data points
    plot(xb, zb, 'k')

    title('Bed profile')
    legend('Sea level', 'Sensors', 'Location', 'SouthEast')
    xlim(xlims)
    ylim([-8, 2])
    xlabel('x [m]')
    ylabel('z [m]')

    % save fig
    print(['figures/62_model', num2str(i)], '-dpng', '-r300')
end


%% ADDITIONAL FIGURE
% figure
% subplot(4,1,1)
% box on
% grid on
% hold on
% for i = 1:n_t
%     plot(xb, Ur(:,i))
%     scatter(xp16, Ur_eg(:,i), '+')
% end

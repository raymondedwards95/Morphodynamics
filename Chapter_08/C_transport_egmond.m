% Sediment transport modelling: Application to the Egmond field work
% Chapter 8.2
% Script 1
%
close all
clear all


%% SETTINGS
D50 = 0.225; % sediment size [mm]
D90 = D50;
rho_s = 2650; % sediment density [kg/m3]
porosity = 0.6; % bed porosity
delta_t = 10 * 24*60*60; % duration [s]
rho_w = 1000; % water density [kg/m3]

% tides
n_t = 3;
tides = {'Low tide', 'Mid tide', 'High tide'};
i_t = 2; % tide to use

% model
% Minimal water depth for computation (we stop the computation when h<hmin)
hmin = 0.2;

% visualizations
xlims = [4200, 5200];
ylims_a = [-3.2, 3.2] .* power(10, -5); % Q
ylims_b = [-0.2, 8]; % O
ylims_c = [-1, 1]; % dz
ylims_d = [-1.5, 1.5] .* power(10, -6); % dQ
ylims_e = [-8, 2]; % zb


%% LOAD DATA
% load H13, T13, Zeta, theta
load('data/52_ParametersEgmond.mat', 'H13', 'T13', 'Zeta', 'theta');

% only use one tide:
H13 = H13(i_t);
T13 = T13(i_t);
Zeta = Zeta(i_t);
theta = theta(i_t);

bedprofile = load('prof1018.txt');
xb = bedprofile(:,1);
zb = bedprofile(:,2);


%% PREPARE
n_c = 5; % number of different cases
n_x = length(xb); % number of points

% convert parameters
Hrms0 = H13 / sqrt(2); % H_{1/3} ~ sqrt{2} H_{rms}
T0 = T13; % T ~ T_{1/3}

% visualizations
if xlims(2) > max(xb) % set xlim, depending on given value and bed profile
    xlims(2) = max(xb);
end


% pre-allocation
Qsx = NaN*ones(n_x, n_c); % save data for all points for all cases
Qsy = Qsx;  Rh = Qsx;
Occ = Qsx;  Oct = Qsx;      Otc = Qsx;  Ott = Qsx;
Uw = Qsx;   Ur = Qsx;       Sk = Qsx;   As = Qsx;
r = Qsx;    phi = Qsx;      R = Qsx;    beta = Qsx;
Urms = Qsx; dQsx_dx = Qsx;  dz = Qsx;   Ux = Qsx;


%% CALCULATIONS
for i = 1:n_c % loop over all different cases
    disp(['Calculating case ', num2str(i), ' of ', num2str(n_c), '. '])

    % set Hrms0 for 8.2.3
    if i == 3
        disp('  Hrms0: half')
        Hrms0__ = Hrms0 ./ 2;
    elseif i == 4
        disp('  Hrms0: double')
        Hrms0__ = Hrms0 .* 2;
    else
        disp('  Hrms0: normal')
        Hrms0__ = Hrms0;
    end

    % set ripples for 8.2.4
    if i == 5
        disp('  Ripples: yes')
        ripples = 1;
    else
        disp('  Ripples: no')
        ripples = 0;
    end

    % run model
    wavedata = BJmodel(Hrms0__, T0, Zeta, theta, bedprofile, hmin);

    % extract variables
    k = wavedata.k;
    h = wavedata.ht;
    Hrms = wavedata.Hrms;

    % find last non NaN
    N_last = find(~isnan(Hrms), 1, 'last');

    % set undertow for 8.2.2 and 8.2.3
    if (2 <= i) && (i <= 4)
        disp('  Undertow: yes')
        % multiply by 100 to get cm/s
        Ux(:,i) = 100 * undertow_magnitude(wavedata.E, wavedata.Er, ...
                                           wavedata.c, wavedata.ht, ...
                                           rho_w, wavedata.Hrms);
    else
        disp('  Undertow: no')
        Ux(:,i) = zeros(n_x, 1);
    end

    % compute variables
    [Uw(:,i)] = velocity_amplitude(h, Hrms, T0);
    [Ur(:,i)] = ursell_number(k, h, Hrms);
    [Sk(:,i), As(:,i)] = sk_as_ursell(Ur(:,i));

    for j = 1:N_last
        [r(j,i)] = computation_r(Sk(j,i), As(j,i));
        [phi(j,i)] = computation_phi(Sk(j,i), As(j,i));
        [u, t] = waveshape(r(j,i), phi(j,i), Uw(j,i), T0);
        [Urms(j,i)] = 100 * rms(u); % multiply by 100 to get cm/s
        [R(j,i), beta(j,i)] = vsk_ask(u, t);

        % sediment transport
        [Qsx(j,i), Qsy(j,i), Occ(j,i), Oct(j,i), Ott(j,i), Otc(j,i), Rh(j,i)] = ...
            SANTOSSmodel(D50, D90, rho_s, T0, Urms(j,i), R(j,i), ...
                         beta(j,i), ripples, Ux(j,i));
    end

    % sediment transport gradient
    dQsx_dx(:,i) = gradient(Qsx(:,i), xb);

    % bed change
    dz(:,i) = -porosity * dQsx_dx(:,i) * delta_t;
end
clear k h Hrms u t Hrms0__ ripples wavedata N_last i j


%% VISUALISATION 1
figure('Name', '821 Sediment transport Egmond')
sgtitle('Sediment transport at Egmond')

% transport
subplot(3,1,1)
hold on; box on; grid on
plot(xb, Qsx(:,1))
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('Cross-shore sediment trasport')
ylabel('Q_{s,x} [m^2/s]')
xlim(xlims)
ylim(ylims_a)
set(gca, 'xticklabel', []) % remove xticklabels

% gradient
subplot(3,1,2)
hold on; box on; grid on
plot(xb, dQsx_dx(:,1))
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('Gradient of cross-shore sediment transport')
ylabel('dQ_{s,x}/dx [m/s]')
xlim(xlims)
ylim(ylims_d)
set(gca, 'xticklabel', []) % remove xticklabels

% bed profile
subplot(3,1,3)
hold on; box on; grid on
plot(xb, zb, 'k')
plot(xlims, Zeta*ones(2,1), '--b') % plot water level

title('Bed profile')
xlabel('x [m]')
ylabel('z [m]')
xlim(xlims)
ylim(ylims_e)
legend('Bed', 'Sea level', 'Location', 'NorthWest')

% save figure
print('figures/82_sediment_transport', '-dpng')


%% VISUALISATION 2
figure('Name', '821 Bed profile changes')
sgtitle('Bed profile change at Egmond')

% change of bed
subplot(2,1,1)
hold on; box on; grid on
plot(xb, dz(:,1))
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('Change of bed elevation')
ylabel('\delta z [m]')
xlim(xlims)
ylim(ylims_c)
set(gca, 'xticklabel', []) % remove xticklabels

% old bed versus new bed
subplot(2,1,2)
hold on; box on; grid on
plot(xb, zb + dz(:,1))
plot(xb, zb, '--k')

title('Bed profile')
xlabel('x [m]')
ylabel('z [m]')
xlim(xlims)
ylim(ylims_e)
legend('New profile', 'Initial profile', 'Location', 'NorthWest')

% save figure
print('figures/82_bedprofile', '-dpng')


%% VISUALISATION 3
figure('Name', '822 Waves and Current')
sgtitle('Effect of undertow')

% Qsx
subplot(3,2,1)
hold on; box on; grid on
for i = 1:2
    plot(xb, Qsx(:,i))
end
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('Q_{s,x}')
ylabel('Q_{s,x} [m^2/s]')
xlim(xlims)
ylim(ylims_a)
% set(gca, 'xticklabel', []) % remove xticklabels
xlabel('x [m]')

% Qsy
subplot(3,2,2)
hold on; box on; grid on
for i = 1:2
    plot(xb, Qsy(:,i))
end
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('Q_{s,y}')
ylabel('Q_{s,y} [m^2/s]')
xlim(xlims)
ylim(ylims_a)
set(gca, 'xticklabel', []) % remove xticklabels
legend('No undertow', 'Undertow')

% Occ
subplot(3,2,3)
hold on; box on; grid on
for i = 1:2
    plot(xb, Occ(:,i))
end
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('\Omega_{c,c}')
ylabel('\Omega_{c,c} [-]')
xlim(xlims)
ylim(ylims_b)
set(gca, 'xticklabel', []) % remove xticklabels

% Oct
subplot(3,2,4)
hold on; box on; grid on
for i = 1:2
    plot(xb, Oct(:,i))
end
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('\Omega_{c,t}')
ylabel('\Omega_{c,t} [-]')
xlim(xlims)
ylim(ylims_b)
set(gca, 'xticklabel', []) % remove xticklabels

% Ott
subplot(3,2,5)
hold on; box on; grid on
for i = 1:2
    plot(xb, Ott(:,i))
end
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('\Omega_{t,t}')
ylabel('\Omega_{t,t} [-]')
xlim(xlims)
ylim(ylims_b)
xlabel('x [m]')

% Otc
subplot(3,2,6)
hold on; box on; grid on
for i = 1:2
    plot(xb, Otc(:,i))
end
plot(xlims, zeros(2,1), '--k') % plot y = 0

title('\Omega_{t,c}')
ylabel('\Omega_{t,c} [-]')
xlim(xlims)
ylim(ylims_b)
xlabel('x [m]')

% save figure
print('figures/82_components_undertow', '-dpng')


%% VISUALISATION 4
for j = 1:4
    if j == 1
        indices = [1];
        name_text = '821 Normal';
        title_text = 'Sediment transport at Egmond';
        labels = {'Waves'};
    elseif j == 2
        indices = [1, 2];
        name_text = '822 Undertow';
        title_text = 'Sediment transport at Egmond (undertow)';
        labels = {'No undertow', 'Undertow'};
    elseif j == 3
        indices = [2, 3, 4];
        name_text = '823 Wave height';
        title_text = 'Sediment transport at Egmond (waveheight)';
        labels = {'Normal', 'Half', 'Double'};
    elseif j == 4
        indices = [1, 5];
        name_text = '824 Ripples';
        title_text = 'Sediment transport at Egmond (ripples)';
        labels = {'No ripples', 'Ripples'};
    else
        indices = [1];
        name_text = '821 Normal';
        title_text = 'Sediment transport at Egmond';
        labels = {'Waves'};
    end

    figure('Name', name_text)
    sgtitle(title_text)

    %% transport
    subplot(3,1,1)
    hold on; box on; grid on
    for i = indices
        plot(xb, Qsx(:,i))
    end
    plot(xlims, zeros(2,1), '--k') % plot y = 0

    title('Cross-shore sediment trasport')
    ylabel('Q_{s,x} [m^2/s]')
    xlim(xlims)
    ylim(ylims_a)
    set(gca, 'xticklabel', []) % remove xticklabels
    legend(labels, 'Location', 'NorthWest')

    %% gradient
    subplot(3,1,2)
    hold on; box on; grid on
    for i = indices
        plot(xb, dQsx_dx(:,i))
    end
    plot(xlims, zeros(2,1), '--k') % plot y = 0

    title('Gradient of cross-shore sediment transport')
    ylabel('dQ_{s,x}/dx [m/s]')
    xlim(xlims)
    ylim(ylims_d)
    set(gca, 'xticklabel', []) % remove xticklabels

    % %% bed change %%% NOTE: if added, change subplot numbers to (4,1,*)
    % subplot(4,1,3)
    % hold on; box on; grid on
    % for i = indices
    %     plot(xb, dz(:,i))
    % end
    % plot(xlims, zeros(2,1), '--k') % plot y = 0

    % title('Change of bed elevation')
    % ylabel('\delta z [m]')
    % xlim(xlims)
    % ylim(ylims_c)
    % set(gca, 'xticklabel', []) % remove xticklabels

    %% bed profile
    subplot(3,1,3)
    hold on; box on; grid on
    plot(xb, zb, '--k')
    plot(xlims, Zeta*ones(2,1), '--b') % plot water level
    for i = indices
        plot(xb, zb + dz(:,i))
    end

    title('Bed profile')
    xlabel('x [m]')
    ylabel('z [m]')
    xlim(xlims)
    ylim(ylims_e)
    legend('Initial bed', 'Sea level', 'Location', 'NorthWest')

    %% save figure
    print(['figures/82_Transport_Egmond_', num2str(j)], '-dpng')
end

clear indices name_text title_text

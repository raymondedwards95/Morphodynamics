% Non-linear wave transformation: Cross-shore evolution of orbital velocity
% Chapter 6.3
% Script 2
%
close all
clear all


%% SETTINGS
% tides
n_t = 3;
tides = {'Low tide', 'Mid tide', 'High tide'};
tide = 1; % make figures for tide 1 (low tide)

% time-series of orbital velocity at following positions
xv = [1000, 4400, 4500, 4700, 4920];

% model
hmin = 0.2; % Minimal water depth for computation (we stop the computation when h<hmin)

% sensors
xp16 = [4478, 4765, 4790, 4814, 4835]; % location of sensors 1,3-6

% visualizations
xlims = [4000, 5200];

%% LOAD DATA
load('data/52_ParametersEgmond.mat', 'H13', 'T13', 'Zeta', 'theta');
load('data/61_WaveShapesEg.mat', 'Sk_eg', 'As_eg');

bedprofile = load('prof1018.txt');
xb = bedprofile(:,1);
zb = bedprofile(:,2);


%% PREPARE
n_x = length(xb); % number of points
n_v = length(xv); % number of locations to determine velocities

% convert parameters
Hrms = H13 / sqrt(2); % H_{1/3} ~ sqrt{2} H_{rms}
T0 = T13; % T ~ T_{1/3}

% pre-allocation
Ur = zeros(n_x, n_t);
Sk = zeros(n_x, n_t);
As = zeros(n_x, n_t);
Uw = zeros(n_x, n_t);

orbit_v = zeros(1000, n_v); % 1000 is the number of elements from function wave_shape
orbit_t = zeros(1000, n_v);


% visualizations
if xlims(2) > max(xb) % set xlim, depending on given value and bed profile
    xlims(2) = max(xb);
end
labels = strings(1, n_v); % used for legends


%% CALCULATIONS
% cross-shore evolution
for i = n_t:-1:1
    wavedata(i) = BJmodel(Hrms(i), T0(i), Zeta(i), theta(i), bedprofile, hmin); % calculate wave evolution

    Ur(:,i) = ursell_number(wavedata(i).k, wavedata(i).ht, wavedata(i).Hrms); % calculate the ursell number for all points

    [Sk(:,i), As(:,i)] = sk_as_ursell(Ur(:,i));

    Uw(:,i) = velocity_amplitude(wavedata(i).ht, wavedata(i).Hrms, T0(i));
end

% time-series
for i = 1:n_v
    Sk_x = interp1(xb, Sk(:,tide), xv(i));
    As_x = interp1(xb, As(:,tide), xv(i));
    Uw_x = interp1(xb, Uw(:,tide), xv(i));

    r_x = computation_r(Sk_x, As_x);
    phi_x = computation_phi(Sk_x, As_x);

    [orbit_v(:,i), orbit_t(:,i)] = waveshape(r_x, phi_x ,Uw_x, T0(1));

    labels(i) = ['x = ', num2str(xv(i))]; % create labels for legends
end
clear Sk_x As_x r_x phi_x Uw_x


%% VISUALIZATION 1
% cross-shore evolution of velocity amplitude
figure('Name', '63 Cross-shore evolution of velocity amplitude')

% velocity amplitude
subplot(2,1,1)
box on
grid on
hold on
% for i = 1:n_v
%     plot(xv(i)*ones(2,1), [-10, 10], '-', 'LineWidth', 0.2) % lines of interest
% end
plot(xb, Uw(:,1))

title('Orbital velocity amplitude')
xlim(xlims)
ylim([0, 1.1])
ylabel('U_w [m/s]')
set(gca, 'xticklabel', []) % remove xticklabels

% bed profile
subplot(2,1,2)
box on
grid on
hold on
% for i = 1:n_v
%     plot(xv(i)*ones(2,1), [-10, 10], '-', 'LineWidth', 0.2) % lines of interest
% end
plot(xlims, Zeta(tide)*ones(2,1), '-.b') % show water level
plot(xb, zb, 'k')

title('Bed profile')
% legend(['Sea level', labels], 'Location', 'SouthEast')
xlim(xlims)
ylim([-10, 2])
xlabel('x [m]')
ylabel('z [m]')

% save figure
print('figures/63_evolution_velocity', '-dpng', '-r300')


%% VISUALIZATION 2
% horizontal orbital velocity at specific points
figure('Name', '63 Horizontal orbital velocities')
box on
grid on
hold on
for i = 1:n_v
    % plot(orbit_t(:,i)/max(orbit_t(:)), orbit_v(:,i)/max(orbit_v(:,1)))
    plot(orbit_t(:,i)/max(orbit_t(:)), orbit_v(:,i))
end
plot([0, 1], [0, 0], '--k') % horizontal line at y=0
% legend(split(num2str(xv), '  '))
legend(labels)
xticks(0:1/8:1) % set xticks
yticks(-1.5:0.25:1.5) % set yticks
ylim([-1.5, 1.5])
xlabel('t/T')
ylabel('u [m/s]')

% save figure
print('figures/63_compare_velocity', '-dpng', '-r300')

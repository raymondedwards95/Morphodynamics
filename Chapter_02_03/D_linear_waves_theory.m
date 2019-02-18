% Linear wave theory: Theoretical study
% Chapter 3.1
% Script 1
%
% Wavelenght and velocity in shallow water and deep water limits
clear all
close all


%% SETTINGS
T = [6, 9, 12]; % list of periods of waves
h = 0:0.5:140; % list of water depths


%% CONSTANTS
g = 9.81; % gravitational acceleration


%% PREPARE CALCULATIONS
n_t = length(T); % number of periods
n_h = length(h); % number of depths

L = zeros(n_h, n_t); % pre-allocate arrays, [rows: depth, coljumns: period]
hL = zeros(n_h, n_t); % pre-allocate arrays
c = zeros(n_h, n_t); % pre-allocate arrays
cg = zeros(n_h, n_t); % pre-allocate arrays


%% CALCULATIONS
for i = 1:n_t % loop over all periods
    for j = 1:n_h % loop over all depths
        L(j,i) = 2 * pi / wave_number(T(i), h(j)); % calculate wavelength
        hL(j,i) = h(j) / L(j,i); % calculate ratio
        c(j,i) = phase_velocity(T(i), h(j)); % calculate phase velocity
        cg(j,i) = group_velocity(T(i), h(j)); % calculate group velocity
    end
end
clear i j


%% FIGURE 1
% figure with wavelengths and ratio h/L
figure
subplot(2,1,1)
hold on
for i = 1:n_t % loop over all periods
    plot(h, L(:,i)) % plot wavelength
end
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s', 'Location', 'SouthEast')
xlabel('h [m]')
ylabel('L [m]')

subplot(2,1,2)
hold on
for i = 1:n_t % loop over all periods
    plot(h, hL(:,i)) % plot ratio h/L
end
plot(h, 0.05*ones(n_h, 1), 'k', 'LineWidth', 1) % plot horizontal lines for shallow water
plot(h, 0.5*ones(n_h, 1), 'k', 'LineWidth', 1) % plot horizontal lines for deep water
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s', 'Location', 'SouthEast')
xlabel('h [m]')
ylabel('h/L')
ylim([0, 0.6])

saveas(gcf, 'figures/3_1_length.png')


%% FIGURE 2
% figure with c and cg
figure
subplot(2,1,1)
hold on
for i = 1:n_t % loop over all periods
    plot(h, c(:,i)) % plot c
end
plot(h, sqrt(g*h)) % plot sqrt(g*h) (c in shallow water limit)
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s', 'Location', 'SouthEast')
xlabel('h [m]')
ylabel('c [m/s]')
ylim([0, 20])

subplot(2,1,2)
hold on
for i = 1:n_t % loop over all periods
    plot(h, cg(:,i)) % plot cg
end
plot(h, sqrt(g*h)) % plot sqrt(g*h) (cg in shallow water limit)
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s')
xlabel('h [m]')
ylabel('c_g [m/s]')
ylim([0, 20])

saveas(gcf, 'figures/3_1_velocity.png')

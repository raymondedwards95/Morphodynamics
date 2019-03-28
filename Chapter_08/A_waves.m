% Sediment transport modelling: Theoretical study
% Preliminary computations
% Chapter 8.1
% Script 1
%
close all
clear all


%% SETTINGS
T = 6; % period (s)
Uw = 1; % velocity amplitude (m/s)
r = [0, 0.6, 0.6, 0.6]; % non-linearity
phi = [0, 0, -pi/2, -pi/4]; % phase


%% PREPARE CALCULATIONS
assert(length(r) == length(phi));
n_c = length(r); % number of combinations of r and phi
n_s = 1000; % number of elements from waveshape
u = zeros(n_s, n_c);
a = u;
t = u;
R = zeros(1, n_c);
Beta = R;


%% CALCULATIONS
for i = 1:n_c % loop over all combinations
    % horizontal orbital velocity
    [u(:,i), t(:,i)] = waveshape(r(i), phi(i), Uw, T);
    % horizontal orbital acceleration
    a(:,i) = acceleration(u(:,i), t(:,i));
    % velocity and acceleration skewness
    [R(i), Beta(i)] = vsk_ask(u(:,i), t(:,i));
end


%% VISUALISATION
for i = 1:n_c
    [t0, t1] = bounds(t(:,i));

    figure('Name', ['81 Velocity ', num2str(i)])
    % sgtitle({['r = ', num2str(r(i)), '; \phi = ', num2str(phi(i))], ...
    %          ['R = ', num2str(R(i)), '; \beta = ', num2str(Beta(i))]})
    sgtitle(['r = ', num2str(r(i)), '; \phi = ', num2str(phi(i))])

    subplot(2,1,1)
    hold on
    box on
    grid on
    plot(t(:,i), u(:,i))
    plot([t0, t1], [0, 0], '--k') % plot horizontal line

    title('Velocity')
    % xlabel('t [s]')
    ylabel('v [m/s]')
    set(gca, 'xticklabel', []) % remove xticklabels
    xlim([t0, t1])

    subplot(2,1,2)
    hold on
    box on
    grid on
    plot(t(:,i), a(:,i))
    plot([t0, t1], [0, 0], '--k') % plot horizontal line

    title('Acceleration')
    xlabel('t [s]')
    ylabel('a [m/s]')
    xlim([t0, t1])

    % save fig
    print(['figures/81_velocity_', num2str(i)], '-dpng')
end

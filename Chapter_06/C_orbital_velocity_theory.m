% Non-linear wave transformation: Cross-shore evolution of orbital velocity
% Chapter 6.3
% Script 1
%
close all
clear all


%% SETTINGS
T = 10; % s
Uw = 1; % m/s

phi_1 = repmat([-pi/2], 1, 3);
phi_2 = repmat([0], 1, 3);
phi_3 = [-pi/2, -pi/4, 0];
phi = cat(1, phi_1, phi_2, phi_3);
clear phi_1 phi_2 phi_3

r_1 = [0, 0.3, 0.6];
r_2 = [0, 0.3, 0.6];
r_3 = repmat([0.6], 1, 3);
r = cat(1, r_1, r_2, r_3);
clear r_1 r_2 r_3

labels = strings(3,3);
labels(1,:) = {'r = 0', 'r = 0.3', 'r = 0.6'};
labels(2,:) = {'r = 0', 'r = 0.3', 'r = 0.6'};
labels(3,:) = {'\Phi = -\pi/2', '\Phi = -\pi/4', '\Phi = 0'};
titles = {'\Phi = -\pi/2', '\Phi = 0', 'r = 0.6'};


%% CALCULATIONS
% calculate waveshape
for j = 1:3 % loop over groups
    for i = 1:3 % loop over cases
        [u(:,j,i), t(:)] = waveshape(r(j,i), phi(j,i), Uw, T);
    end
end

% scale (normalize)
u_ = u/Uw;
t_ = t/T;


%% VISUALIZATIONS
for j = 1:3 % loop over groups
    figure('Name', ['63 Orbital velocity theory ', num2str(j)])
    box on
    grid on
    hold on
    for i = 1:3 % loop over cases
        plot(t_, u_(:,j,i))
    end
    plot([0, 1], [0, 0], '--k') % horizontal line at y=0
    xlabel('t/T')
    ylabel('u/U_w')
    xticks(0:1/8:1)
    ylim([-1.5, 1.5])
    title(['Orbital velocity'; titles(j)])
    legend(labels(j,:))

    % save figure
    print(['figures/63_theory_velocity_', num2str(j)], '-dpng', '-r300')
end

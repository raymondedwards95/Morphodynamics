% Sediment transport modelling: Theoretical study
% Sediment transport analysis
% Chapter 8.1
% Script 2
%
close all
clear all


%% SETTINGS
T = 7; % period [s]
Uw = 1.2; % orbital velocity amplitude [m/s]
r = 0:0.15:0.6; % non-linearity parameter
phi = [-pi/2, 0]; % phase [-]
D50 = [0.3, 0.1]; % sediment size [mm]
D90 = D50;
Rhos = 2650; % sediment density [kg/m3]
ripples = 0; % include ripples? [-]
Ux = 0; % undertow magnitude [m/s]


%% PREPARE
n_r = length(r);
n_p = length(phi);
n_d = length(D50);

% pre-allocation
labels = strings(n_p, n_d);

Qsx = zeros(n_r, n_p, n_d);
Qsy = Qsx;
Occ = Qsx; Oct = Qsx;
Ott = Qsx; Otc = Qsx;

u = zeros(1000, n_p, n_r); % waveshape returns 1000 elements for u
t = u;
R = zeros(n_p, n_r);
beta = R;
Urms = R;

% visualisations
[r_min, r_max] = bounds(r);

%% COMPUTE
for i = 1:n_r % loop over all r
    for j = 1:n_p % loop over all phi
        % 1: computation of the time-series of orbital velocity
        [u(:,j,i), t(:,j,i)] = waveshape(r(i), phi(j), Uw, T);

        % 2: computation of the velocity skewness R and
        %    the acceleration skewness beta
        [R(j,i), beta(j,i)] = vsk_ask(u(:,j,i), t(:,j,i));

        % 3: computation of the root-mean squared orbital velocity Urms:
        %    should be in CM/S!
        Urms(j,i) = rms(u(:,j,i)) * 100; 
        % Urms(j,i) = sqrt(mean(power(u(:,j,i), 2))) * 100;
        % Urms(j,i) = std(u(:,j,i)) * 100;

        % 4: sediment transport calculation
        for k = 1:n_d % loop over all sediment sizes
            [Qsx(i,j,k), Qsy(i,j,k), Occ(i,j,k), Oct(i,j,k), Ott(i,j,k), Otc(i,j,k)] = SANTOSSmodel(D50(k), D90(k), Rhos, T, Urms(j,i), R(j,i), beta(j,i), ripples, Ux);
        end
    end
end


%% VISUALISATION 1
figure('Name', '81 Waves and transport')

% Qsx
subplot(3,2,1)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Qsx(:,j,k), '-o', 'MarkerSize', 5)
        labels(k,j) = ['D50 = ', num2str(D50(k)), '; \phi = ', num2str(phi(j))];
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

legend(labels(:), 'Location', 'SouthEast')
xlabel('r')
title('Q_{s,x}')
ylabel('Q_{s,x}')


% Qsy
subplot(3,2,2)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Qsy(:,j,k), '-o', 'MarkerSize', 5)
        labels(k,j) = ['D50 = ', num2str(D50(k)), '; \phi = ', num2str(phi(j))];
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

legend(labels(:), 'Location', 'SouthEast')
xlabel('r')
title('Q_{s,y}')
ylabel('Q_{s,y}')


% Occ
subplot(3,2,3)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Occ(:,j,k), '-o', 'MarkerSize', 5)
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

xlabel('r')
title('\Omega_{c,c}')
ylabel('\Omega_{c,c}')


% Oct
subplot(3,2,4)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Oct(:,j,k), '-o', 'MarkerSize', 5)
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

xlabel('r')
title('\Omega_{c,t}')
ylabel('\Omega_{c,t}')


% Ott
subplot(3,2,5)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Ott(:,j,k), '-o', 'MarkerSize', 5)
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

xlabel('r')
title('\Omega_{t,t}')
ylabel('\Omega_{t,t}')


% Otc
subplot(3,2,6)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Otc(:,j,k), '-o', 'MarkerSize', 5)
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

xlabel('r')
title('\Omega_{t,c}')
ylabel('\Omega_{t,c}')


% save figure
print('figures/81_transport', '-dpng')


%% VISUALISATION 2
figure('Name', '81 Waves and transport (2)')

% main plot
subplot(3,2,1)
hold on; grid on; box on
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        plot(r, Qsx(:,j,k), '-o', 'MarkerSize', 5)
        labels(k,j) = ['D50 = ', num2str(D50(k)), '; \phi = ', num2str(phi(j))];
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

legend(labels(:), 'Location', 'best')
xlabel('r')
title('Q_{s,x}')
ylabel('Q_{s,x}')

% main plot 2
subplot(3,2,2)
hold on; grid on; box on
for j = 1:n_p
    for k = 1:n_d
        plot(r, Qsy(:,j,k), '-o', 'MarkerSize', 5)
        labels(k,j) = ['D50 = ', num2str(D50(k)), '; \phi = ', num2str(phi(j))];
    end
end
plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

legend(labels(:), 'Location', 'best')
xlabel('r')
title('Q_{s,y}')
ylabel('Q_{s,y}')

% sub plot
for j = 1:n_p % loop over all phi
    for k = 1:n_d % loop over all D50
        subplot(3,2,2+2*(j-1)+k)
        hold on; grid on; box on
        plot(r, Occ(:,j,k), '-o', 'MarkerSize', 5)
        plot(r, Oct(:,j,k), '-o', 'MarkerSize', 5)
        plot(r, Ott(:,j,k), '-o', 'MarkerSize', 5)
        plot(r, Otc(:,j,k), '-o', 'MarkerSize', 5)
        plot([r_min, r_max], [0, 0], '--k') % plot horizontal line

        xlabel('r')
        title(labels(k,j))
        ylabel('\Omega')
        legend('\Omega_{c,c}', '\Omega_{c,t}', '\Omega_{t,t}', '\Omega_{t,c}', 'Location', 'best')
    end
end

% save figure
print('figures/81_transport_alt', '-dpng')

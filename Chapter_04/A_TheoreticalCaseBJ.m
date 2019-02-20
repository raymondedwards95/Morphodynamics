% Main script for the computation of the cross-shore wave evolution 
% for a simplified bar-trough case based on the Battjes and Janssen (1978) model 

%------------------------------------
%           Initialisation
%------------------------------------

clear all;
close all;

% Definition of cross-shore coordinates (m)
x = (1:1:500)';  

% Definition of zb, bed elevation relative to mean water level (m)
zb = -9 * ones(size(x));   
zb = zb + (x<=260).*(x-160)/20  + ((x<=300)&(x>260)).*(-1/40*(x-260) + 5) + (x>300).*(1/20*(x-300) + 4); 

% Definition of the array profile, input argument for BJmodel
profile = [x zb];

% Offshore wave conditions
Hrms0 = 1;      % Root mean square wave height (m)
theta0 = 0;     % Angle of incidence (degrees)
T0 = 10;        % Characteristic period (s)
Zeta = 0;       % Mean water level (m)

% Model parameter 
hmin = 0.2;     % Minimal water depth for computation
                % (we stop the computation when h<hmin)

%------------------------------------------------
%       Computation of wave characteristics
%------------------------------------------------

waves = BJmodel(Hrms0, T0, Zeta, theta0, profile, hmin);


%------------------------------------
%           Visualisation 
%------------------------------------
x_lims = [0, 450];

figure
subplot(4,1,1)
plot(waves.x, waves.Hrms)
ylabel('Hrms [m]')
xlim(x_lims)
ylim([0, 2.5])
title('Root mean square wave height')

subplot(4,1,2)
hold on
plot(waves.x, waves.eta)
plot(waves.x, waves.eta(1)*ones(length(waves.x)), '--k')
hold off
ylabel('SSE [m]')
xlim(x_lims)
ylim([-0.2, 0.2])
title('Sea surface elevation')

subplot(4,1,3)
hold on
plot(waves.x, waves.Dbr)
plot(waves.x, waves.Dr)
hold off
legend('D_{br}', 'D{r}', 'Location', 'NorthWest')
ylabel('D [W/m^2]')
xlim(x_lims)
ylim([0, 550])
title('Dissipation')

subplot(4,1,4)
hold on
plot(waves.x, waves.z, 'k')
plot(waves.x, Zeta*ones(size(x)), '-.')
hold off
xlabel('x [m]')
ylabel('zb [m]')
xlim(x_lims)
ylim([-20, 5])
title('Bed profile')

saveas(gcf, ['figures/4_2_theory_', num2str(Zeta), '_', num2str(theta0), '_', num2str(Hrms0), '.png'])

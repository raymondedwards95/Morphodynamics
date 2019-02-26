% Modelling cross-shore wave transformation: Theoretical example
% Chapter 4.2
% script 1
%
clear all
close all


%% SETTINGS and PARAMETERS
% offshore wave conditions
% NOTE: script will do all combinations of parameters
% T0: number: Characteristic period (s)
T0 = 10; % default: 10;
% Hrms0: list: Root mean square wave height (m)
Hrms0 = [1]; % default: [1];   alternative: [0.5, 1, 2];
% theta0: list: Angle of incidence (degrees)
theta0 = [0]; % default: [0];   alternative: [0, 22.5, 45];
% Zeta: list: Mean water level (m)
Zeta = [-1, 0, 1]; % default: [0];   alternative: [-1, 0, 1];

% Model parameter 
hmin = 0.2; % Minimal water depth for computation (we stop the computation when h<hmin)

% Definition of cross-shore coordinates (m)
x = (1:1:500)';  

% Definition of zb, bed elevation relative to mean water level (m)
zb = -9 * ones(size(x));   
zb = zb + (x<=260).*(x-160)/20  + ((x<=300)&(x>260)).*(-1/40*(x-260) + 5) + (x>300).*(1/20*(x-300) + 4); 

% Definition of the array profile, input argument for BJmodel
profile = [x zb];


%% PREPARE CALCULATIONS
n_Hrms = length(Hrms0); % number of parameters
n_theta = length(theta0); % number of parameters
n_Zeta = length(Zeta); % number of parameters

labels = strings(n_Hrms, n_theta, n_Zeta); % pre-allocate labels for legends in figures


%% COMPUTATIONS
% NOTE: loops are backwards for pre-allocation purposes
for h = n_Hrms:-1:1 % loop over all given initial Hrms
    for t = n_theta:-1:1 % loop over all given initial theta
        for z = n_Zeta:-1:1% loop over all given Zeta
            disp([num2str(h), ' ', num2str(t), ' ', num2str(z)])
            data(h,t,z) = BJmodel(Hrms0(h),T0,Zeta(z),theta0(t),profile,hmin);
            labels(h,t,z) = ['\zeta = ', num2str(Zeta(z)), ' m; H_{rms,0} = ', num2str(Hrms0(h)), ' m; \theta_0 = ', num2str(theta0(t)), ''];
        end
    end
end


%% FIGURES
x_lims = [0, 450];

figure

% Hrms
subplot(5,1,1)
hold on
y_lim2 = 0;
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).Hrms)
            y_lim2 = max([y_lim2, max(data(h,t,z).Hrms)]); % find max, for plot limits
        end
    end
end
hold off
box on
ylabel('Hrms [m]')
xlim(x_lims)
ylim([0, 1.1*y_lim2])
title('Root mean square wave height')
set(gca, 'xticklabel', []) % remove xticklabels

% eta
subplot(5,1,2)
hold on
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).eta)
        end
    end
end
plot(data(1).x, zeros(length(data(1).x)), '--k') % horizontal line at eta=0
hold off
box on
ylabel('SSE [m]')
xlim(x_lims)
ylim([-0.05, 0.05])
title('Sea surface elevation')
set(gca, 'xticklabel', [])

% Dbr
subplot(5,1,3)
hold on
y_lim2 = 0;
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).Dbr)
            y_lim2 = max([y_lim2, max(data(h,t,z).Dbr)]); % find max, for plot limits
        end
    end
end
hold off
box on
ylabel('D_{br} [W/m^2]')
xlim(x_lims)
ylim([0, 1.1*y_lim2])
title('Dissipation: breaking')
set(gca, 'xticklabel', [])

% Dr
subplot(5,1,4)
hold on
y_lim2 = 0;
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).Dr)
            y_lim2 = max([y_lim2, max(data(h,t,z).Dr)]); % find max, for plot limits
        end
    end
end
hold off
box on
ylabel('D_{r} [W/m^2]')
xlim(x_lims)
ylim([0, 1.1*y_lim2])
title('Dissipation: roller')
set(gca, 'xticklabel', [])
legend(labels, 'Location', 'West')

% bed profile
subplot(5,1,5)
hold on
plot(x, zb, 'k') % bed
plot(x, Zeta.*ones(size(x)), '-.') % sea level
hold off
box on
xlabel('x [m]')
ylabel('zb [m]')
xlim(x_lims)
% ylim([-20, 5])
title('Bed profile')

set(gcf,'position',[1, 1, 800, 2000]) % x0, y0, width, height (pixels)

saveas(gcf, ['figures/4_2_theory_compare_', num2str(n_Zeta), '_', num2str(n_theta), '_', num2str(n_Hrms), '.png'])

clear h t z

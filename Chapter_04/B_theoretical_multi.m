%
% Chapter 4.1
% script 2
%
clear all
close all


%% SETTINGS and PARAMETERS
% offshore wave conditions
T0 = 10;            % number: Characteristic period (s)
Hrms0 = [0.5, 1, 2, 3];        % list: Root mean square wave height (m)
theta0 = [0];       % list: Angle of incidence (degrees)
Zeta = [0];  % list: Mean water level (m)

% Model parameter 
hmin = 0.2;         % Minimal water depth for computation (we stop the computation when h<hmin)

% Definition of cross-shore coordinates (m)
x = (1:1:500)';  

% Definition of zb, bed elevation relative to mean water level (m)
zb = -9 * ones(size(x));   
zb = zb + (x<=260).*(x-160)/20  + ((x<=300)&(x>260)).*(-1/40*(x-260) + 5) + (x>300).*(1/20*(x-300) + 4); 

% Definition of the array profile, input argument for BJmodel
profile = [x zb];


%% PREPARE CALCULATIONS
n_Hrms = length(Hrms0);
n_theta = length(theta0);
n_Zeta = length(Zeta);

labels = strings(n_Hrms, n_theta, n_Zeta); % labels for legends in figures


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
subplot(5,1,1)
hold on
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).Hrms)
        end
    end
end
hold off
ylabel('Hrms [m]')
xlim(x_lims)
% ylim([0, 2.5])
title('Root mean square wave height')
set(gca, 'xticklabel', [])

subplot(5,1,2)
hold on
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).eta)
        end
    end
end
plot(data(1).x, zeros(length(data(1).x)), '--k')
hold off
ylabel('SSE [m]')
xlim(x_lims)
% ylim([-0.2, 0.2])
title('Sea surface elevation')
set(gca, 'xticklabel', [])

subplot(5,1,3)
hold on
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).Dbr)
        end
    end
end
hold off
ylabel('D_{br} [W/m^2]')
xlim(x_lims)
% ylim([0, 550])
title('Dissipation - breaking')
set(gca, 'xticklabel', [])

subplot(5,1,4)
hold on
for h = 1:n_Hrms
    for t = 1:n_theta
        for z = 1:n_Zeta
            plot(data(h,t,z).x, data(h,t,z).Dr)
        end
    end
end
hold off
ylabel('D_{r} [W/m^2]')
xlim(x_lims)
% ylim([0, 550])
title('Dissipation - rolling')
set(gca, 'xticklabel', [])
legend(labels, 'Location', 'West')

subplot(5,1,5)
hold on
plot(x, zb, 'k')
plot(x, Zeta.*ones(size(x)), '-.')
hold off
xlabel('x [m]')
ylabel('zb [m]')
xlim(x_lims)
% ylim([-20, 5])
title('Bed profile')

set(gcf,'position',[1, 1, 800, 2000]) % x0, y0, width, height (pixels)

saveas(gcf, ['figures/4_2_theory_compare_', num2str(n_Zeta), '_', num2str(n_theta), '_', num2str(n_Hrms), '.png'])


% Wave-by-wave analysis of the low tide data, Egmond Coast 3d dataset
% The objective of the script is to compute wave statistics at different
% locations along a cross-shore transect and examine their cross-shore
% evolution


% -------------------------------------
%            Initialisation
% -------------------------------------
clear all
close all

% load data
data = load('lowTide.txt');

% constants
g = 9.81;                    % acceleration of gravity (m/s^2)
rho = 1025;                  % water density (kg/m^3)
fs = 2;                      % sampling frequency (Hz)
Npos = 5;                    % number of cross-shore positions considered

% Initialisation vectors 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(Npos,1);  % root mean square height (m)
H13_tot  = zeros(Npos,1);  % significant wave height (m)
Hm_tot   = zeros(Npos,1);  % mean wave height (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------

for i=1:Npos  % loop on the positions
    
    wavedata = zero_crossing(data(:,i), fs);
    Hrms_tot(i,1) = rms_height(wavedata(:,1));  
    H13_tot(i,1) = significant_height(wavedata(:,1));
    Hm_tot(i,1) = mean(wavedata(:,1));
    
end

% --------------------------------------
%            Positions and bed
% --------------------------------------

pos_x = [4478 4765 4790 4814 4835];
pos_name = ['P1' 'P3' 'P4' 'P5' 'P6'];

bedprofile = load('prof1018.txt');
bed_x = bedprofile(:,1);
bed_y = bedprofile(:,2);

% --------------------------------------
%                  Output
% --------------------------------------

x_left = 4400;
x_right = max(pos_x)+10;

% visualisation of outputs
figure;
subplot(211)
hold on
plot(pos_x, Hm_tot, 'blue') % lines
plot(pos_x, H13_tot, 'green')
plot(pos_x, Hrms_tot, 'red')
scatter(pos_x, Hm_tot, 'blue') % dots
scatter(pos_x, H13_tot, 'green')
scatter(pos_x, Hrms_tot, 'red')
hold off
legend('H_{m}', 'H_{1/3}', 'H_{rms}')
xlim([x_left, x_right])
ylim([0, 2])
title('Statistics of low-tide')
xlabel('Distance from bouy (m)')
ylabel('Height (m)')

subplot(212)
hold on
plot(bed_x, bed_y) % line
scatter(pos_x, interp1(bed_x, bed_y, pos_x)) % dots, using interp to get y-value for x of all positions
hold off
xlim([x_left, x_right])
ylim([-7, 0])
title('Bed profile')
xlabel('Distance from bouy (m)')
ylabel('Height (m)')

% save data
saveas(gcf, 'figures/1_2_lowtides.png')

save('lowtidestatisticsdata', 'Hm_tot', 'Hrms_tot', 'H13_tot')
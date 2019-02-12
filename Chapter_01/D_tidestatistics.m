% Wave-by-wave analysis of the low tide data, Egmond Coast 3d dataset
% The objective of the script is to compute wave statistics at different
% locations along a cross-shore transect and examine their cross-shore
% evolution


% -------------------------------------
%            Initialisation
% -------------------------------------
clear all
close all

% Load data
data01 = load('lowTide.txt');
data02 = load('midTide.txt');
data03 = load('highTide.txt');
data=[data01,data02,data03];

bedprofile = load('prof1018.txt');

% constants
g = 9.81;                    % Acceleration of gravity (m/s^2)
rho = 1025;                  % Water density (kg/m^3)
fs = 2;                      % Sampling frequency (Hz)
Npos = 5;                    % Number of cross-shore positions considered
Ndata = 3;                   % Number of tide types (low, mid, high)

% Initialisation vectors 
% These vectors will be used to store the wave statistics at each position
Hrms_tot = zeros(Npos,Ndata);  % Root mean square height, H_rms (m)
H13_tot  = zeros(Npos,Ndata);  % Significant wave height, H_1/3 (m)
Hm_tot   = zeros(Npos,Ndata);  % Mean wave height, H_m (m)

% --------------------------------------
%     Computation of wave statistics
% --------------------------------------
for j=1:Ndata % Loop on the type of tide
   for i=1:Npos  % Loop on the pressure sensor positions
      dummylist = zero_crossing(data(:,i+Npos*(j-1)),fs);
      Hm_tot(i,j) = mean(dummylist(:,1));
      H13_tot(i,j) = significant_height(dummylist(:,1));
      Hrms_tot(i,j) = rms_height(dummylist(:,1));
   end
end

% --------------------------------------
%                  Output
% --------------------------------------

% Fitting the data
x = reshape(Hrms_tot,1,15);     % All H_rms values (m)
y = reshape(H13_tot,1,15);      % All H_13 values (m)
z = reshape(Hm_tot,1,15);       % All h_mean values (m)
fit = polyfit(x,y,1);           % Linear fit of H_rms as function of H_1/3
fit2 = polyfit(x,z,1);          % Linear fit of H_m as function of H_rms
plotfit = polyval(fit,x);       

% Visualisation of outputs
figure, scatter(Hrms_tot(:,1),H13_tot(:,1))
hold on, scatter(Hrms_tot(:,2),H13_tot(:,2))
hold on, scatter(Hrms_tot(:,3),H13_tot(:,3))
hold on, plot(x,plotfit)
legend('Low tide','Mid tide','High tide','Linear fit','location','southeast')
title('Height comparison')
xlabel('H_{rms} (m)')
ylabel('H_{1/3} (m)')

% Saving the data
save('StatisticsEgmond','Hm_tot','Hrms_tot','H13_tot');

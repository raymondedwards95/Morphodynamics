% Sediment transport observations: Analysis of the hydrodynamics
% Chapter 7.1
% Script 1
%
close all
clear all


%% SETTINGS
n_b = 8; % divide the timeseries by n_b parts
n_s = 4096; % number of samples

f_1 = 0.005; % lowest frequency infragravity
f_2 = 0.05; % division infragravity - sea swell


%% LOAD DATA
load('ucData.mat', 'condition1', 'condition2', 'condition3');
conditions = [condition1, condition2, condition3];
clear condition1 condition2 condition3


%% PREPARATIONS
nfft = round(n_s / n_b);

t = repmat(transpose(0:n_s-1), 1, 3); % time

n_c = length(conditions);

eta = zeros(n_s, n_c);
eta_low = eta;
eta_high = eta;

f_s = zeros(1, n_c);
f_n = f_s;
h = f_s;

Hm0 = zeros(1, n_c);
Hm0_ss = Hm0;
Hm0_inf = Hm0;

Hm0_ss_h = Hm0;
Hm0_inf_Hm0_ss = Hm0;

S = zeros(round(nfft/2+1), n_c);
f = S;

Sk = zeros(1, n_c);
As = Sk;
Sk_high = Sk;
As_high = Sk;

u = eta;
u_low = eta;
u_high = eta;
u_mean = Sk;
u_mean_sum = eta;

c = eta;
c_low = eta;
c_high = eta;
c_mean = u_mean;
c_mean_sum = u_mean_sum;

q = u_mean;
q_low = u_mean;
q_high = u_mean;
q_mean = u_mean;
q_mean_sum = u_mean;


%% RESHAPE and RESAMPLE DATA
for i = 1:n_c
    eta(:,i) = conditions(i).data(:,1);
    u(:,i) = conditions(i).data(:,2);
    c(:,i) = conditions(i).data(:,3);

    f_s(1,i) = conditions(i).Fs;
    f_n(1,i) = f_s(i) / 2;
    h(1,i) = conditions(i).h;

    % find time for timeseries
    t(:,i) = t(:,i) / f_s(i);

    % resample eta
    eta_low(:,i) = fft_filter(eta(:,i), f_s(i), f_1, f_2); % remove high freqs
    eta_high(:,i) = fft_filter(eta(:,i), f_s(i), f_2, f_n(i)); % remove low freqs

    % resample u
    u_low(:,i) = fft_filter(u(:,i), f_s(i), f_1, f_2); % remove high freqs
    u_high(:,i) = fft_filter(u(:,i), f_s(i), f_2, f_n(i)); % remove low freqs

    % resample c
    c_low(:,i) = fft_filter(c(:,i), f_s(i), f_1, f_2); % remove high freqs
    c_high(:,i) = fft_filter(c(:,i), f_s(i), f_2, f_n(i)); % remove low freqs
end
clear conditions


%% COMPUTATIONS
for i = 1:n_c
    % spectrum
    [S(:,i), f(:,i), ~, ~] = VarianceDensitySpectrum(eta(:,i), nfft, f_s(1,i));

    % wave heights (from moment) % H = 4 sqrt(m)
    % inf: f_1 = 0.005 < f < 0.05 = f_2 Hz
    % ss:  f_2 = 0.05  < f < 1.   = f_n  Hz
    Hm0(i) = 4 * sqrt(spectral_moment(f(:,i), S(:,i), f_1, f_n(i), 0));
    Hm0_inf(i) = 4 * sqrt(spectral_moment(f(:,i), S(:,i), f_1, f_2, 0));
    Hm0_ss(i) = 4 * sqrt(spectral_moment(f(:,i), S(:,i), f_2, f_n(i), 0));

    % calculate ratios
    Hm0_ss_h(i) = Hm0_ss(i) / h(i);
    Hm0_inf_Hm0_ss(i) = Hm0_inf(i) / Hm0_ss(i);

    % skewness and asymmetry
    [Sk(i), As(i)] = sk_as(eta(:,i));
    [Sk_high(i), As_high(i)] = sk_as(eta_high(:,i));

    % velocity
    u_mean(i) = mean(u(:,i));
    u_mean_sum(:,i) = u_mean(i) + u_high(:,i) + u_low(:,i);

    % concentration
    c_mean(i) = mean(c(:,i));
    c_mean_sum(:,i) = c_mean(i) + c_high(:,i) + c_low(:,i);

    % transport
    q(i) = mean(u(:,i) .* c(:,i)); % = mean(u c)
    q_mean(i) = u_mean(i) * c_mean(i); % = mean(u) mean(c)
    q_high(i) = mean(u_high(:,i) .* c_high(:,i)); % mean(u_h c_h)
    q_low(i) = mean(u_low(:,i) .* c_low(:,i)); % mean(u_l c_l)
    q_mean_sum(i) = q_mean(i) + q_high(i) + q_low(i);
end


%% VISUALIZATION spectrum
figure('Name', ['71 Spectrum', num2str(i)])
for i = 1:n_c
    subplot(n_c,1,i)
    plot(f(:,i), S(:,i))
    grid on
    box on

    xlabel('f [Hz]')
    ylabel('S [m^2/Hz]')
end


%% VISUALIZATION velocity
figure('Name', '71 Velocity decomposition check')
plot(t, u_mean_sum(:,1) - u(:,1))
plot([min(t(:,i)), max(t(:,i))], [0, 0], '--k') % horizontal line at y=0


%% VISUALIZATION velocity: high freq vs low freq
figure('Name', '71 Velocity: low vs high frequency')
sgtitle('Cross-shore velocity')

for i = 1:n_c
    subplot(n_c,1,i)
    box on
    grid on
    hold on
    plot(t(:,i), u_high(:,i))
    plot(t(:,i), u_low(:,i))
    plot([min(t(:,i)), max(t(:,i))], [0, 0], '--k') % horizontal line at y=0

    title(['Condition ', num2str(i)])
    xlim([min(t(:,i)), max(t(:,i))])
    ylim([-1.5, 1.5])
    xlabel('t [s]')
    ylabel('u [m/s]')
end
legend('u_{hf}', 'u_{lf}')


%% VISUALIZATION sediment
for i = 1:n_c
    figure('Name', ['72 Preliminary analysis ', num2str(i)])
    sgtitle(['Condition ', num2str(i)])

    % velocity
    subplot(2,1,1)
    box on
    grid on
    hold on
    plot(t(:,i), u(:,i)-u_mean(i))
    plot(t(:,i), u_low(:,i))
    plot([min(t(:,i)), max(t(:,i))], [0, 0], '--k') % horizontal line at y=0
    
    title('Cross-shore velocity')
    legend('u - mean(u)', 'u_{lf}')
    xlim([min(t(:,i)), max(t(:,i))])
    ylabel('u [m/s]')
    set(gca, 'xticklabel', []) % remove xticklabels

    % concentration
    subplot(2,1,2)
    box on
    grid on
    hold on
    plot(t(:,i), c(:,i))

    title('Concentration')
    xlim([min(t(:,i)), max(t(:,i))])  
    xlabel('t [s]')
    ylabel('c [kg/m^3]')  
end

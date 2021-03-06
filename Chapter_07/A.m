% Sediment transport observations
% Chapter 7
% Script 1
%
close all
clear all


%% SETTINGS
n_b = 8; % divide the timeseries by n_b parts, number of blocks used in spectrum is equal to 2*n_b-1, for n_b = 8 -> 15 blocks

f_1 = 0.005; % lowest frequency infragravity
f_2 = 0.05; % separation infragravity - sea swell


%% LOAD DATA
load('ucData.mat', 'condition1', 'condition2', 'condition3');
conditions = [condition1, condition2, condition3]; % reshape data, to be able to use loops
clear condition1 condition2 condition3

titles = {'Condition 1', 'Condition 2', 'Condition 3'};


%% PREPARATIONS
n_c = length(conditions); % number of conditions
n_s = length(conditions(1).data); % total number of samples, assuming that all time-series have the same number

nfft = round(n_s / n_b); % number of samples per block

t = repmat(transpose(0:n_s-1), 1, 3); % time-axis

% pre-allocation
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

S = zeros(round(nfft/2+1), n_c); % spectrum for each condition
f = S;

Sk = zeros(1, n_c);
Sk_low = Sk;
Sk_high = Sk;
As = Sk;
As_low = Sk;
As_high = Sk;

u = zeros(n_s, n_c);
u_low = u;
u_high = u;
u_sum = u;
u_mean = zeros(1, n_c);
u_error = u_mean;
u_diff = u_mean;

c = u;
c_low = u;
c_high = u;
c_sum = u_sum;
c_mean = u_mean;
c_error = u_mean;
c_diff = u_mean;
corr_uc = zeros(3, n_c);

q = u_mean;
q_low = u_mean;
q_high = u_mean;
q_mean = u_mean;
q_sum = u_mean;
q_error = u_mean;
q_diff = u_mean;


%% RESHAPE and FILTER DATA
for i = 1:n_c % loop over all conditions
    % save all data in separate variables
    % time-series
    eta(:,i) = conditions(i).data(:,1); % sea-surface elevation
    u(:,i) = conditions(i).data(:,2); % velocity
    c(:,i) = conditions(i).data(:,3); % concentration

    % numbers
    f_s(i) = conditions(i).Fs; % sampling frequency
    f_n(i) = f_s(i) / 2; % Nyquist frequency
    h(i) = conditions(i).h; % water depth

    % find (exact) time for timeseries
    t(:,i) = t(:,i) / f_s(i);

    % filter eta
    eta_low(:,i) = fft_filter(eta(:,i), f_s(i), f_1, f_2); % remove high freqs
    eta_high(:,i) = fft_filter(eta(:,i), f_s(i), f_2, f_n(i)); % remove low freqs

    % filter u
    u_low(:,i) = fft_filter(u(:,i), f_s(i), f_1, f_2); % remove high freqs
    u_high(:,i) = fft_filter(u(:,i), f_s(i), f_2, f_n(i)); % remove low freqs

    % filter c
    c_low(:,i) = fft_filter(c(:,i), f_s(i), f_1, f_2); % remove high freqs
    c_high(:,i) = fft_filter(c(:,i), f_s(i), f_2, f_n(i)); % remove low freqs
end
clear conditions


%% COMPUTATIONS
for i = 1:n_c % loop over all conditions
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
    [Sk_low(i), As_low(i)] = sk_as(eta_low(:,i));
    [Sk_high(i), As_high(i)] = sk_as(eta_high(:,i));

    % velocity
    u_mean(i) = mean(u(:,i));
    u_sum(:,i) = u_mean(i) + u_high(:,i) + u_low(:,i);
    u_error(i) = rms_error(u_sum(:,i), u(:,i));
    u_diff(i) = mean(abs(u_sum(:,i) - u(:,i)));

    % concentration
    c_mean(i) = mean(c(:,i));
    c_sum(:,i) = c_mean(i) + c_high(:,i) + c_low(:,i);
    c_error(i) = rms_error(c_sum(:,i), c(:,i));
    c_diff(i) = mean(abs(c_sum(:,i) - c(:,i)));
    corr_1 = corrcoef(c(:,i), u(:,i));
    corr_2 = corrcoef(c(:,i), u_low(:,i));
    corr_3 = corrcoef(c(:,i), u_high(:,i));
    corr_uc(:,i) = [corr_1(1,2), corr_2(1,2), corr_3(1,2)];

    % transport
    q(i) = mean(u(:,i) .* c(:,i)); % = mean(u c)
    q_mean(i) = u_mean(i) * c_mean(i); % = mean(u) mean(c)
    q_high(i) = mean(u_high(:,i) .* c_high(:,i)); % mean(u_h c_h)
    q_low(i) = mean(u_low(:,i) .* c_low(:,i)); % mean(u_l c_l)
    q_sum(i) = q_mean(i) + q_high(i) + q_low(i);
    q_error(i) = rms_error(q_sum(i), q(i));
    q_diff(i) = mean(abs(q_sum(i) - q(i)));
end
clear corr_1 corr_2 corr_3


%% VISUALIZATION spectrum
figure('Name', '71 Spectrum')
for i = 1:n_c % loop over all conditions
    subplot(n_c,1,i)
    % yyaxis left
    hold on
    box on
    grid on
    plot(f(:,i), S(:,i))

    % yyaxis right
    % semilogy(f(:,i), S(:,i))

    title(['Condition ', num2str(i)])
    xlabel('f [Hz]')
    ylabel('S [m^2/Hz]')
    xlim([0, 0.5])
    ylim([0, 1.05*max(S(:,i))])
    if max(S(:,i)) < 0.5
        yticks(0:0.1:max(S(:,i)))
    elseif max(S(:,i)) < 1.5
        yticks(0:0.2:max(S(:,i)))
    elseif max(S(:,i)) < 5
        yticks(0:1:max(S(:,i)))
    elseif max(S(:,i)) < 10
        yticks(0:2:max(S(:,i)))
    end
end

% save figure
print('figures/71_spectrum', '-dpng')


%% VISUALIZATION velocity decomposition
figure('Name', '71 Velocity decomposition check')
hold on
box on
grid on
for i = 1:n_c
    plot(t(:,i), (u_sum(:,i) - u(:,i)) ./ max(u(:,i)))
end
plot([min(t(:,i)), max(t(:,i))], [0, 0], '--k') % horizontal line at y=0

title('Relative difference between decomposed and original u')
xlabel('t [s]')
ylabel('Relative difference')
xlim([min(t(:,i)), max(t(:,i))])
legend(titles)

% save figure
print('figures/71_differences', '-dpng')


%% VISUALIZATION velocity: high freq vs low freq
figure('Name', '71 Velocity: low vs high frequency')
sgtitle('Cross-shore velocity')

for i = 1:n_c % loop over conditions
    subplot(n_c,1,i)
    hold on
    box on
    grid on
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

% save figure
print('figures/71_decompose_velocity', '-dpng')


%% VISUALIZATION sediment
for i = 1:n_c % loop over all conditions
    figure('Name', ['72 Preliminary analysis ', num2str(i)])
    sgtitle(['Condition ', num2str(i)])

    % velocity
    subplot(2,1,1)
    hold on
    box on
    grid on
    plot(t(:,i), u(:,i)-u_mean(i))
    plot(t(:,i), u_low(:,i))
    plot([min(t(:,i)), max(t(:,i))], [0, 0], '--k') % horizontal line at y=0

    title('Cross-shore velocity')
    legend('u - u_{mean}', 'u_{lf}')
    xlim([min(t(:,i)), max(t(:,i))])
    ylabel('u [m/s]')
    set(gca, 'xticklabel', []) % remove xticklabels

    % concentration
    subplot(2,1,2)
    hold on
    box on
    grid on
    plot(t(:,i), c(:,i))

    title('Concentration')
    xlim([min(t(:,i)), max(t(:,i))])
    xlabel('t [s]')
    ylabel('c [kg/m^3]')

    % save figure
    print(['figures/72_velocity_concentration_', num2str(i)], '-dpng')
end


%% VISUALIZATION sediment
figure('Name', '72 Transport values')
hold on
box on
grid on
bar(categorical(titles), transpose([q; q_mean; q_low; q_high; q_sum]), 'grouped')

title('Sediment transport')
ylabel('q [kg/m^2/s]')
legend('Total', 'Mean', 'Low frequency', 'High frequency', 'Sum')

% save figure
print('figures/72_transport', '-dpng')

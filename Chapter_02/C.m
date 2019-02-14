% Spectral analysis: Computation of spectral wave characteristics
% Chapter 2.2
% Script 1
%
% Comparing block sizes using one location during low tide
clear all
close all

f_s = 2;
f_N = f_s / 2;

fmin = [0, 0.005, 0.05];
fmax = [f_N, 0.05, f_N];
flabels = {'m_{0}'; 'm_{0,inf}'; 'm_{0,ss}'};

blocks = 15;

lowtide = load('lowTide.txt');
midtide = load('midTide.txt');
hightide = load('highTide.txt');

data = cat(3, lowtide, midtide, hightide);
clear lowtide midtide hightide

bedprofile = load('prof1018.txt');
x_b = bedprofile(:,1);
y_b = bedprofile(:,2);
x_p = [4478, 4765, 4790, 4814, 4835];

[n_s, n_p, n_t] = size(data);
assert(length(fmin) == length(fmax), 'fmin and fmax have different lengths');
n_f = min([length(fmin), length(fmax)]);

nfft = round(n_s ./ ((blocks + 1) / 2));

% Hm0_all = zeros(n_p, n_t, n_f);
m0_all = zeros(n_p, n_t, n_f);
H13 = zeros(n_p, n_t);
fp = zeros(n_p, n_t);

for i = 1:n_p % loop over all positions
    for j = 1:n_t % loop over all tides
        temp_data = data(:,i,j);
        [wave] = zero_crossing(temp_data, f_s);
        [S, f, e, c] = VarianceDensitySpectrum(temp_data, nfft, f_s);
        for k = 1:n_f % loop over all pairs of fmin and fmax
            m0_all(i,j,k) = spectral_moment(f, S, fmin(k), fmax(k), 0);
        end
        H13(i,j) = significant_height(wave);

        [m, mi] = max(S);
        fp(i,j) = f(mi);
    end
end
clear m mi S f e c temp_data i j k wave

Tp = 1 ./ fp;
Hm0_all = 4 * sqrt(m0_all);

Hm0 = Hm0_all(:,:,1);
Hm0_ss = Hm0_all(:,:,3);
Hm0_inf = Hm0_all(:,:,2);

%%% FIGURE 1
Hm0f = reshape(Hm0, n_t*n_p, 1);
H13f = reshape(H13, n_t*n_p, 1);
coeff = polyfit(Hm0f, H13f, 1);
figure;
hold on
scatter(Hm0f, H13f)
plot(Hm0f, polyval(coeff, Hm0f))
hold off
xlabel('H_{m0}')
ylabel('H_{1/3}')
legend('data', ['fit: slope = ', num2str(coeff(1))], 'Location', 'NorthWest')
saveas(gcf, 'figures/2_2_spectral_compare.png')

%%% FIGURE 2
tide = 1;
ylims = [2, 0.5, 2];
x_left = 4350;
x_right = max(x_p)+10;

figure;
sgtitle('Spectral wave heights')

for k = 1:n_f
    subplot(n_f+1, 1, k);
    hold on
    scatter(x_p, Hm0_all(:,tide,k), '+')
    plot(x_p, Hm0_all(:,tide,k), '--')
    hold off
    xlim([x_left, x_right])
    ylim([0, ylims(k)])
    set(gca, 'xticklabel', [])
    ylabel([flabels(k), ' [m]'])
    legend([flabels(k)], 'Location', 'NorthWest')
end

subplot(n_f+1, 1, n_f+1)
hold on
plot(x_b, y_b) % line
scatter(x_p, interp1(x_b, y_b, x_p)) % dots, using interp to get y-value for x of all positions
hold off
xlim([x_left, x_right])
ylim([-7, 0])
title('Bed profile')
xlabel('Distance from bouy [m]')
ylabel('Height [m]')
saveas(gcf, 'figures/2_2_spectral_low.png')

save('data_spectral', 'Hm0', 'Hm0_ss', 'Hm0_inf', 'Tp', 'H13')
clear k ylims x_left x_right
% Spectral analysis: Computation of density spectrum, part 1
% Chapter 2.1
% Script 1
%
% Comparing block sizes using one location during low tide
clear all
close all


%%% SETTINGS
f_s = 2; % sampling frequency
blocks = [1, 3, 7, 15, 31]; % specify total number of blocks, including overlapping blocks


%%% READ DATA
lowtide = load('lowTide.txt'); % data


%%% PREPARE CALCULATIONS
[n_s, n_p] = size(lowtide); % number of samples and number of sensor positions
n_b = length(blocks);
nfft = round(n_s ./ ((blocks + 1) / 2)); % length of each block

S = zeros(n_s, n_b); % save S for plots, for each block size, there is a list (1-dim array)
f = zeros(n_s, n_b); % save f for plots
edf = zeros(1, n_b);
conf95Interval = zeros(2, n_b); % save for plots, 2 numbers for each block size


%%% CALCULATIONS
for i=1:1:n_b % loop over all block sizes
    % NOTE: complex indices for S and f due to fact that not all results are of the same length
    % Length of S_i is always (nfft/2 + 1)
    [S(1:nfft(i)/2+1,i), f(1:nfft(i)/2+1,i), edf(i), conf95Interval(:,i)] = VarianceDensitySpectrum(lowtide(:,1), nfft(i), f_s);
end
clear i


%%% FIGURE
figure;
sgtitle('Variance density spectrum for P1')
for i=1:1:n_b % loop over all block lengths
    subplot(n_b, 1, i)
    hold on
    plot(f(:,i), S(:,i)) % plot center
    plot(f(:,i), S(:,i) * conf95Interval(1,i), '-r', 'LineWidth', 0.1) % plot confidence interval
    plot(f(:,i), S(:,i) * conf95Interval(2,i), '-r', 'LineWidth', 0.1)
    hold off
    if i == n_b % check if subplot is last subplot
        xlabel('Frequency [Hz]')
    else % subplot is not last subplot
        set(gca, 'xticklabel', [])
    end
    xlim([0, 0.7])
    ylabel('S [m^2 Hz^{-1}]')
    legend([num2str(blocks(i)), ' blocks'])
end
clear i
saveas(gcf, 'figures/2_1_p1.png')

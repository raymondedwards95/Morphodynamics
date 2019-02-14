% Spectral analysis: Computation of density spectrum, part 2
% Chapter 2.1
% Script 2
%
% Comparing different locations during low tide
% Most code (adapted) from chapter 2.1 script 1
clear all
close all


%%% SETTINGS
f_s = 2;
blocks = [15];


%%% READ DATA
lowtide = load('lowTide.txt');
names = {'P1'; 'P3'; 'P4'; 'P5'; 'P6'};


%%% PREPARE CALCULATIONS
[n_s, n_p] = size(lowtide);
n_b = length(blocks);
nfft = round(n_s ./ ((blocks + 1) / 2));

S = zeros(n_s, n_b, n_p);
f = zeros(n_s, n_b, n_p);
edf = zeros(1, n_b, n_p);
conf95Interval = zeros(2, n_b, n_p);


%%% CALCULATIONS
for j=1:1:n_p % loop over all positions
    for i=1:1:n_b % loop over all block sizes
        % NOTE: complex indices for S and f due to fact that not all results are of the same length
        [S(1:nfft(i)/2+1,i,j), f(1:nfft(i)/2+1,i,j), edf(i,j), conf95Interval(:,i,j)] = VarianceDensitySpectrum(lowtide(:,j), nfft(i), f_s);
    end
end


%%% FIGURE
for i=1:1:n_b % loop over all blocks
    figure;
    sgtitle(['Variance density spectrum with ', num2str(blocks(i)), ' blocks']) % sgtitle --> title for figure, rather than for subplot
    for j=1:1:n_p % loop over all locations
        subplot(n_p,1,j)
        hold on
        plot(f(:,i,j), S(:,i,j)) % plot center
        % plot(f(:,i,j), S(:,i,j) * conf95Interval(1,i,j)) % plot confidence interval
        % plot(f(:,i,j), S(:,i,j) * conf95Interval(2,i,j))
        hold off
        if j == n_p % check if subplot is last subplot
            xlabel('Frequency [Hz]')
        else % subplot is not last subplot
            set(gca, 'xticklabel', [])
        end
        ylabel('S [m^2 Hz^{-1}]')
        legend(names(j))
        xlim([0, 0.7])
        ylim([0, 2])
    end
    saveas(gcf, ['figures/2_1_all_', num2str(blocks(i)),'.png'])
end

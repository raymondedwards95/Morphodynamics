clear all
close all

f_s = 2;

lowtide = load('lowTide.txt');

[n_s, m_s] = size(lowtide);

blocks = [1, 3, 7, 15, 31];
n_b = length(blocks);
nfft = round(n_s ./ ((blocks + 1) / 2));

S = zeros(n_s, n_b);
f = zeros(n_s, n_b);
edf = zeros(1, n_b);
conf95Interval = zeros(2, n_b);

for i=1:1:n_b
    [S(1:nfft(i)/2+1,i), f(1:nfft(i)/2+1,i), edf(i), conf95Interval(:,i)] = VarianceDensitySpectrum(lowtide(:,1),nfft(i),f_s);
end

figure;
sgtitle('Variance density spectrum for P1')
for i=1:1:n_b
    subplot(n_b,1,i)
    hold on
    plot(f(:,i), S(:,i))
    plot(f(:,i), S(:,i) * conf95Interval(1,i), '-r', 'LineWidth', 0.1)
    plot(f(:,i), S(:,i) * conf95Interval(2,i), '-r', 'LineWidth', 0.1)
    hold off
    if i == n_b
        xlabel('Frequency [Hz]')
    else
        set(gca, 'xticklabel', [])
    end
    xlim([0, 0.7])
    ylabel('S [m^2 Hz^{-1}]')
    legend([num2str(blocks(i)), ' blocks'])
end
saveas(gcf, 'figures/2_1_p1.png')

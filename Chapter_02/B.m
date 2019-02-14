clear all
close all

f_s = 2;

lowtide = load('lowTide.txt');

[n_s, n_p] = size(lowtide);
names = ['P1'; 'P3'; 'P4'; 'P5'; 'P6'];

blocks = [15];
n_b = length(blocks);
nfft = round(n_s ./ ((blocks + 1) / 2));

S = zeros(n_s, n_b, n_p);
f = zeros(n_s, n_b, n_p);
edf = zeros(1, n_b, n_p);
conf95Interval = zeros(2, n_b, n_p);

for j=1:1:n_p
    for i=1:1:n_b
        [S(1:nfft(i)/2+1,i,j), f(1:nfft(i)/2+1,i,j), edf(i,j), conf95Interval(:,i,j)] = VarianceDensitySpectrum(lowtide(:,j),nfft(i),f_s);
    end
end

for i=1:1:n_b
    figure;
    sgtitle(['Variance density spectrum with ', num2str(blocks(i)), ' blocks'])
    for j=1:1:n_p
        subplot(n_p,1,j)
        hold on
        plot(f(:,i,j), S(:,i,j))
        % plot(f(:,i,j), S(:,i,j) * conf95Interval(1,i,j))
        % plot(f(:,i,j), S(:,i,j) * conf95Interval(2,i,j))
        hold off
        if j == n_p
            xlabel('Frequency [Hz]')
        else
            set(gca, 'xticklabel', [])
        end
        ylabel('S [m^2 Hz^{-1}]')
        legend(['Location ', names(j,:)])
        xlim([0, 0.7])
        ylim([0, 2])
    end
    saveas(gcf, ['figures/2_1_all_', num2str(blocks(i)),'.png'])
end

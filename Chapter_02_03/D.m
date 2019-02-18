% Linear wave theory: Theoretical study
% Chapter 3.1
% Script 1
%
% 
clear all
close all


%%% SETTINGS
T = [6, 9, 12];
h = 0:0.5:140;


%%% CONSTANTS
g = 9.81;


%%% PREPARE CALCULATIONS
n_t = length(T);
n_h = length(h);

L = zeros(n_h, n_t);
hL = zeros(n_h, n_t);
c = zeros(n_h, n_t);
cg = zeros(n_h, n_t);


%%% CALCULATIONS
for i = 1:n_t
    for j = 1:n_h
        L(j,i) = 2 * pi / wave_number(T(i), h(j));
        hL(j,i) = h(j) / L(j,i);
        c(j,i) = phase_velocity(T(i), h(j));
        cg(j,i) = group_velocity(T(i), h(j));
    end
end
clear i j


%%% FIGURE 1
figure
subplot(2,1,1)
hold on
for i = 1:n_t
    plot(h, L(:,i))
end
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s')
xlabel('h [m]')
ylabel('L [m]')

subplot(2,1,2)
hold on
for i = 1:n_t
    plot(h, hL(:,i))
end
plot(h, 0.05*ones(n_h, 1), 'k', 'LineWidth', 1)
plot(h, 0.5*ones(n_h, 1), 'k', 'LineWidth', 1)
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s')
xlabel('h [m]')
ylabel('h/L')
ylim([0, 1])


%%% FIGURE 2
figure
subplot(2,1,1)
hold on
for i = 1:n_t
    plot(h, c(:,i))
end
plot(h, sqrt(g*h))
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s')
xlabel('h [m]')
ylabel('c [m/s]')

subplot(2,1,2)
hold on
for i = 1:n_t
    plot(h, cg(:,i))
end
plot(h, sqrt(g*h))
hold off
legend('T = 6 s', 'T = 9 s', 'T = 12 s')
xlabel('h [m]')
ylabel('c_g [m/s]')

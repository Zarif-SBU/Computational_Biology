n0 = 1;
t2 = 90; %min

%unrestricted growth model
t = 0:180:72*60;
unrestricted = n0 * 2.^(t/t2);
plot(t/60, unrestricted, '-o');
grid on
hold on

% Logistic equation

K = 2*10^11; % Define carrying capacity
n0 = 1; % Initial population
r0 = 0.0077; % Growth rate
t = 0:180:72*60; % Time in minutes

logistic = n0*exp(r0*t).*(K./(K-n0+n0*exp(r0*t)));
plot(t/60, logistic, '-o', 'DisplayName', 'Logistic');

xlabel('Time (hours)');
ylabel('Population');
title('Exponential vs Logistic');
legend('Exponential, N_0 = 1, t2 = 90(min)', 'Logistic, N_0 = 1, R_0 = 0.0077, K = 2*10^{11}')
axis([0,100,0,3.5e14]) %changes range of y axis

hold off
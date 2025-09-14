% Logistic equation
K = 2*10^11; % Define carrying capacity
n0 = 10000; % Initial population
r0 = 0.0077; % Growth rate
t = 0:90:100*60; % Time in minutes

logistic = n0*exp(r0*t).*(K./(K-n0+n0*exp(r0*t)));

plot(t/60, logistic, '-o', 'DisplayName', 'Logistic');
xlabel('Time (1.5 hour increment)');
ylabel('Population');
title('Logistic');
legend('Logistic, N_0 = 10000, R_0 = 0.0077, K = 2*10^{11}', 'Location', 'northeast')
axis([0,100,10000,2.3e11])
grid on


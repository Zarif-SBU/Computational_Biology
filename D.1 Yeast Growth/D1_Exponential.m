%unrestricted growth model
n0 = 10000;
t2 = 90; %min
t = 0:30:45*60;

unrestricted = n0 * 2.^(t/t2);

plot(t/60, unrestricted, '-o');
grid on
xlabel('Time (hours)');
ylabel('Population');
title('Exponential');
legend('Exponential, N_0 = 10000, t2 = 90(min)', 'Location','northwest')
grid on

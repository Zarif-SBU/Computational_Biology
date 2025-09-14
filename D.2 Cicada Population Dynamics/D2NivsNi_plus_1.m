N0 = 100; % Starting population
R0 = 24; % Per capita growth rate
K = 8.75e10; % Carrying capacity
a = R0 / K; % Competitive rate
b = 1;
cycleLen = 17; % Cicada life cycle (in years)
years = 170; % Total years
numCycles = floor(years / cycleLen);

pop = zeros(1, numCycles + 1);
pop(1) = N0; % Init pop

for i = 2:numCycles + 1
    Ni = pop(i - 1);
    pop(i) = (R0 * Ni) / ((1 + a * Ni)^b);
end
Ni = pop(1:end-1);
Ni_plus_1 = pop(2:end);
plot(Ni, Ni_plus_1, '-bo');
hold on
plot(Ni,Ni)
xlabel('N_i');
ylabel('N_i+1');
title('N_i+1 vs N_i');
legend('N_i+1 vs N_i', 'y = x','Location','northwest')
hold off
grid on
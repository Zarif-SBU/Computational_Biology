N0 = 500000 ; % Starting population
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

plot(0:numCycles, pop, '-bo');
xlabel('Cycle Number');
ylabel('Population');
title('Cicada Population Growth Over Time');
legend('Hassel, N0=500000, b=1, a=R0/k, k=8.75e10, R0=24','Location','northwest')
grid on;
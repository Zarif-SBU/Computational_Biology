N = 10000000;
I0 = 4 / N;
R0 = 0;
S0 = 1 - (I0 + R0);

tau = 10; % in days
r0 = 3;
beta = r0 / tau;
gamma = 1 / tau;

total_time = 160; % Total simulation time in days

S = zeros(1, total_time);
I = zeros(1, total_time);
R = zeros(1, total_time);
S(1) = S0;
I(1) = I0;
R(1) = R0;

highest_slope = 0;
highest_increase_day = 0;

for i = 1:total_time
    S(i + 1) = S(i) - beta * S(i) * I(i);
    I(i + 1) = I(i) + (beta * S(i) * I(i) - gamma * I(i));
    R(i + 1) = R(i) + gamma * I(i);

    if (I(i + 1) - I(i)) > highest_slope
        highest_slope = (I(i + 1) - I(i));
        highest_increase_day = (i + 1);
    end
end

% Plotting results
plot(0:1:total_time, S, '-g');
hold on;
plot(0:1:total_time, I, '-r');
plot(0:1:total_time, R, '-b');
xlabel('Days');
ylabel('Fractional Population');
title('COVID-19 SIR Model Simulation');
legend('Suscpetible','Infected', 'Recovered');
grid on;
hold off;

disp(highest_increase_day);
disp(highest_slope);
disp(R(end) + I(end));
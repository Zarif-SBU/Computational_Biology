N = 10000000;
I0 = 4 / N;
R0 = .75;
S0 = 1 - (I0 + R0);

tau = 10; % in days
r0 = 3;
beta = r0 / tau;
gamma = 1 / tau;

dt = 0.1; % Time step size to help make simulation more accurate
total_time = 1000; % Total simulation time in days
time_steps = total_time / dt;

S = zeros(1, time_steps);
I = zeros(1, time_steps);
R = zeros(1, time_steps);
S(1) = S0;
I(1) = I0;
R(1) = R0;

highest_slope = 0;
highest_increase_day = 0;

for i = 1:time_steps
    S(i + 1) = S(i) - beta * S(i) * I(i) * dt;
    I(i + 1) = I(i) + (beta * S(i) * I(i) - gamma * I(i)) * dt;
    R(i + 1) = R(i) + gamma * I(i) * dt;

    if (I(i + 1) - I(i)) / dt > highest_slope
        highest_slope = (I(i + 1) - I(i)) / dt;
        highest_increase_day = (i + 1) * dt;
    end
end

% Plotting results
%plot(0:dt:total_time, S, '-g');
hold on;
plot(0:dt:total_time, I, '-r');
plot(0:dt:total_time, R, '-b');
xlabel('Days');
ylabel('Fractional Population');
title('COVID-19 SIR Model Simulation');
legend('Infected', 'Recovered');
grid on;
hold off;

disp(highest_increase_day);
disp(highest_slope);
disp(R(end) + I(end));
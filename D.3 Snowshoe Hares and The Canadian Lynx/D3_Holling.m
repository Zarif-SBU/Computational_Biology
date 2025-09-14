%D3
global alpha gamma epsilon beta kappa;
alpha = .2205; %
beta = 4.2801; %
gamma = .365; %
epsilon = .015; %Efficiency
kappa = (1/6)/365; %handling time per prey caught

function dU = dU_dt(U, V)
    global alpha gamma kappa
    dU = alpha * U - ((gamma * U) / (1 + gamma * kappa * U)) * V;
end

function dV = dV_dt(U, V)
    global epsilon gamma beta kappa
    dV = epsilon * ((gamma * U) / (1 + gamma * kappa * U)) * V - beta * V;
end

U0 = 400; %init hare pop
V0 = 1; %init lynx pop
h = .001; %Step Size
Ut = zeros(1, 40001);
Vt = zeros(1, 40001);
t_values = zeros(1, 40001);
Ut(1) = U0;
Vt(1) = V0;
t = 0;

for i = 1:40000
    t = t + h;
    Ut(i+1) = Ut(i) + dU_dt(Ut(i), Vt(i)) * h;
    Vt(i+1) = Vt(i) + dV_dt(Ut(i), Vt(i)) * h;
    t_values(i+1) = t; %we i + 1 to skip first index since we know that will 0

end
hold on
plot(t_values, Ut, 'DisplayName', 'Hare');
plot(t_values, Vt, 'DisplayName', 'Lynx');

xlabel('Time (hours)');
ylabel('Population per {km^2}');
title('Hare and Lynx growth over time');
legend('Hare', 'Lynx');

hold off

%discuss your observations, repeat sym with 800, 2 and 200, .5

% Calculate stationary points and nullclines
U_cline = beta / ((epsilon - beta*kappa) * gamma); % Lynx nullcline (U*)
V_cline = (alpha / gamma) + (alpha * kappa * U_cline); % Hare nullcline (V*)

% Phase plane plot (Hare vs Lynx)
figure;
plot(Ut, Vt, 'LineWidth', 1.5); % Hare vs Lynx trajectory
xlabel('Hare Population (per km^2)');
ylabel('Lynx Population (per km^2)');
title('Phase Plane of Hare vs Lynx Populations');
grid on;

% Plot the stationary points
hold on;
plot(0, 0, 'ro', 'MarkerSize', 10, 'DisplayName', 'Stationary Point (0,0)'); % Trivial stationary point
plot(U_cline, V_cline, 'go', 'MarkerSize', 10, 'DisplayName', ['Stationary Point (', num2str(U_cline), ', ', num2str(V_cline), ')']); 

% Plot nullclines
yline(V_cline, 'b--', 'DisplayName', 'Hare Nullcline V* = \alpha/\gamma', 'LineWidth', 1.5);
xline(U_cline, 'm--', 'DisplayName', 'Lynx Nullcline U* = \beta/(\epsilon\gamma)', 'LineWidth', 1.5);

legend('Trajectory', 'Stationary Point (0,0)', ['Stationary Point (', num2str(U_cline), ', ', num2str(V_cline), ')']);
hold off;
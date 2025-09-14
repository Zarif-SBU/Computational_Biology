%D3
global alpha gamma epsilon beta;
alpha = 0.2205; % per capita annual growth rate 
beta = 4.2801; % per capita annual deatrh rate
gamma = .365; % per capita predation rate
epsilon = .015; % Efficiency

function dU = dU_dt(U, V) %Hare equation
    global alpha gamma
    dU = (alpha * U) - (gamma * U * V);         
end

function dV = dV_dt(U, V) %Lynx equation
    global epsilon gamma beta
    dV = (epsilon * gamma * U * V) - (beta * V);
end

U0 = 400; %init hare pop
V0 = 1; %init lynx pop
h = .001; %Step Size
total_time = 40; %hrs
num_steps = total_time/h; %number of steps
Ut = zeros(1, num_steps+1); %initializing array to store our hare population at each step
Vt = zeros(1, num_steps+1); %initializing array to store our lynx population at each step
t_values = zeros(1, num_steps+1); %initializing array to store each step
Ut(1) = U0;
Vt(1) = V0;
t = 0;

max_dU = 0; % maximum rate of change for hare 
max_dV = 0; % maxim rate of change for lynx 
ut = 0
vt = 0
% Euler method for population simulation
for i = 1:num_steps
    t = t + h;
    
    dU = dU_dt(Ut(i), Vt(i));
    dV = dV_dt(Ut(i), Vt(i));
    
    Ut(i+1) = Ut(i) + dU * h;
    Vt(i+1) = Vt(i) + dV * h;
    
    t_values(i+1) = t;
    
    % update the rate of change
    if abs(dU) > max_dU
        max_dU = abs(dU);
    end
    if abs(dV) > max_dV
        max_dV = abs(dV);
    end
end

disp(max_dU)
disp(max_dV)

%code below plots population vs time
hold on
plot(t_values, Ut, 'DisplayName', 'Hare');
plot(t_values, Vt, 'DisplayName', 'Lynx');

xlabel('Time (hours)');
ylabel('Population per {km^2}');
title('Hare and Lynx growth over time');
legend(['U0 = ' num2str(U0) ''], ['V0 = ' num2str(V0)]);

hold off

% Calculate stationary points and nullclines
U_cline = beta / (epsilon * gamma); % Lynx nullcline (vertical)
V_cline = alpha / gamma;            % Hare nullcline (horizontal)

% Phase plane plot
figure;
plot(Ut, Vt, 'LineWidth', 1.5); % Hare vs Lynx trajectory
xlabel('Hare Population (per km^2)');
ylabel('Lynx Population (per km^2)');
title('Phase Plane of Hare vs Lynx Populations');
grid on;

% plots the stationary points
hold on;
plot(0, 0, 'ro', 'MarkerSize', 10, 'DisplayName', 'Stationary Point (0,0)'); % Trivial stationary point
plot(U_cline, V_cline, 'go', 'MarkerSize', 10, 'DisplayName', ['Stationary Point (', num2str(U_cline), ', ', num2str(V_cline), ')']); 

% plots nullclines
yline(V_cline, 'b--', 'DisplayName', 'Hare Nullcline V* = \alpha/\gamma', 'LineWidth', 1.5);
xline(U_cline, 'm--', 'DisplayName', 'Lynx Nullcline U* = \beta/(\epsilon\gamma)', 'LineWidth', 1.5);

legend('Trajectory', 'Stationary Point (0,0)', ['Stationary Point (', num2str(U_cline), ', ', num2str(V_cline), ')']);
hold off;

%discuss your observations, repeat sym with 800, 2 and 200, .5
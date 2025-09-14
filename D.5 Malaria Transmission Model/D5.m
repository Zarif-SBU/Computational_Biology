global beta epsilon1 epsilon2 gamma alpha N_H N_M
beta = 5;
epsilon1 = .01;
epsilon2 = .1;
gamma = .05;
alpha = .02;

N_H = 10^6;
N_M = 10^4;
I_H0 = N_H * 0;
S_H0 = N_H - I_H0;
I_M0 = N_M * 0.1;
S_M0 = N_M - I_M0;

total_time = 200; % days
h = .01; % time step
num_steps = total_time/h;

function dS_H = H_dS_dt(S_H, I_M , I_H)
    global beta epsilon1 gamma N_H
    dS_H = -epsilon1 * beta * (S_H / N_H) * I_M + gamma * I_H;
end

function dI_H = H_dI_dt(S_H, I_M, I_H)
    global beta epsilon1 gamma N_H
    dI_H = epsilon1 * beta * (S_H / N_H) * I_M - gamma * I_H;
end

function dS_M = M_dS_dt(I_H, S_M)
    global alpha N_M epsilon2 beta N_H
    dS_M = alpha * N_M - epsilon2 * beta * (I_H / N_H) * S_M - alpha * S_M;
end

function dI_M = M_dI_dt(I_H, S_M, I_M)
    global alpha N_H epsilon2 beta
    dI_M = epsilon2 *  beta * (I_H / N_H) * S_M - alpha * I_M;
end

H_I = zeros(1, 1 + num_steps);
H_S = zeros(1, 1 + num_steps);
M_I = zeros(1, 1 + num_steps);
M_S = zeros(1, 1 + num_steps);
t_values = zeros(1, num_steps+1); %initializing array to store each step
H_I(1) = I_H0;
H_S(1) = S_H0;
M_I(1) = I_M0;
M_S(1) = S_M0;
t = 0;
for i = 1:num_steps
    t = t + h;
    
    dS_H = H_dS_dt(H_S(i), M_I(i), H_I(i));
    H_S(i+1) = H_S(i) + dS_H * h;

    dI_H = H_dI_dt(H_S(i), M_I(i), H_I(i));
    H_I(i+1) = H_I(i) + dI_H * h;

    dS_M = M_dS_dt(H_I(i), M_S(i));
    M_S(i+1) = M_S(i) + dS_M * h;

    dI_M = M_dI_dt(H_I(i), M_S(i), M_I(i));
    M_I(i+1) = M_I(i) + dI_M * h;

    t_values(i+1) = t;
end

plot(t_values, H_S, 'b-', 'DisplayName', 'Susceptible Humans')
hold on
plot(t_values, H_I, 'r-', 'DisplayName', 'Infected Humans')
plot(t_values, M_S, 'g-', 'DisplayName', 'Susceptible Mosquitoes')
plot(t_values, M_I, 'k-', 'DisplayName', 'Infected Mosquitoes')
xlabel('Time (days)')
ylabel('Population')
title('Malaria Model Simulation')
legend('show')
grid on
hold off
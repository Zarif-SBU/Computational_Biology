% Model Parameters
G_Na = 400; %nS
G_K = 200; %nS
G_L = 2; %nS
E_Na = 99; %mV
E_K = -85; %mV
V_L = -65; %mV
C = 2; %pF

dt = .01; %
total_time = 200; %ms
t = 0:dt:total_time;
num_steps = length(t);
Ie_vals = 101:1:140;
firing_rates = zeros(1, 40);

t0 = 40; % ms 
I0 = 200; % pA
I_e = zeros(1, num_steps);
I_e(t >= t0) = I0;
I_Na = zeros(1, num_steps);
I_K = zeros(1, num_steps);
I_L = zeros(1, num_steps);

V = zeros(1, num_steps);
m = zeros(1, num_steps);
h = zeros(1, num_steps);
n = zeros(1, num_steps);

spike_threshold = -20; % mV
time_between_spikes = 5; % ms
total_spikes = 0;
last_spike = 0;
all_spike_times = [];
all_spike_voltages = [];

%Initialization
V(1) = V_L;

alpha_m = 0.1 * (V(1) + 40) / (1 - exp(-0.1 * (V(1) + 40)));
beta_m = 4 * exp(-0.0556 * (V(1) + 65));
m(1) = alpha_m / (alpha_m + beta_m);

alpha_h = 0.07 * exp(-0.05 * (V(1) + 65));
beta_h = 1 / (1 + exp(-0.1 * (V(1) + 35)));
h(1) = alpha_h / (alpha_h + beta_h);

alpha_n = 0.01 * (V(1) + 55) / (1 - exp(-0.1 * (V(1) + 55)));
beta_n = 0.125 * exp(-0.0125 * (V(1) + 65));
n(1) = alpha_n / (alpha_n + beta_n);

tic
%Eular model
for i = 1:num_steps-1
    alpha_m = .1 * (V(i) + 40)/ (1 - exp(-.1 * (V(i) + 40)));
    beta_m = 4 * exp(-.0556 * (V(i) + 65));

    alpha_h = 0.07 * exp(-.05 * (V(i) + 65));
    beta_h = 1 / (1 + exp(-.1 * (V(i) + 35)));
    
    alpha_n =  0.01 * (V(i) + 55) / (1 - exp(-0.1 * (V(i) + 55)));
    beta_n = .125 * exp(-.0125 * (V(i) + 65));

    dm_dt = alpha_m * (1-m(i)) - beta_m * m(i);
    dh_dt = alpha_h * (1-h(i)) - beta_h * h(i);
    dn_dt = alpha_n * (1-n(i)) - beta_n * n(i);
    
    m(i+1) = m(i) + dm_dt * dt;
    h(i+1) = h(i) + dh_dt * dt;
    n(i+1) = n(i) + dn_dt * dt;

    I_Na(i) = G_Na * m(i)^3 * h(i) * (V(i) - E_Na);
    I_K(i) = G_K * n(i)^4 * (V(i) - E_K);
    I_L(i) = G_L * (V(i) - V_L);

    dV_dt = -I_Na(i) - I_K(i) - I_L(i) + I_e(i);
    V(i+1) = V(i) + (dV_dt/C) * dt;
    if V(i+1) >= spike_threshold && t(i+1) - last_spike > time_between_spikes
        total_spikes = total_spikes + 1;
        all_spike_times = [all_spike_times, t(i+1)];
        all_spike_voltages = [all_spike_voltages, V(i+1)];
        last_spike = t(i+1);
    end
end
toc
I_m = I_Na + I_K + I_L;


figure('Position', [100 100 600 800]);
% Membrane Potential
subplot(5,1,1);
hold on;
plot(t, V);
plot(all_spike_times, all_spike_voltages, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
xlabel('Time (ms)');
ylabel('V (mV)');
title('Membrane Potential');
grid on;
ylim([-100 100]); 
% xlim([72.2,82.2])

% Total Membrane Current 
subplot(5,1,2);
plot(t, I_m);
xlabel('Time (ms)');
ylabel('I_m (pA)');
title('Total Membrane Current');
grid on;
ylim([min(I_m)-50 max(I_m)+50]); 
% xlim([72.2,82.2])

% Sodium Current
subplot(5,1,3);
plot(t, I_Na);
xlabel('Time (ms)');
ylabel('I_{Na} (pA)');
title('Sodium Current');
grid on;
ylim([min(I_Na)-500 max(I_Na)+500]);
% xlim([72.2,82.2])

% Potassium Current
subplot(5,1,4);
plot(t, I_K);
xlabel('Time (ms)');
ylabel('I_K (pA)');
title('Potassium Current');
grid on;
ylim([min(I_K)-500 max(I_K)+500]);
% xlim([72.2,82.2])

% External Current (I_e)
subplot(5,1,5);
plot(t, I_e);
xlabel('Time (ms)');
ylabel('I_e (pA)');
title('External Current');
grid on;
ylim([-50 250]);
% xlim([72.2,82.2])

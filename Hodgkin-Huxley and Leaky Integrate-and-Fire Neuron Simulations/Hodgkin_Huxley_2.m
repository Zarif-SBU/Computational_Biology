% Model Parameters
G_Na = 400; %nS
G_K = 200; %nS
G_L = 2; %nS
E_Na = 99; %mV
E_K = -85; %mV
V_L = -65; %mV
C = 2; %pF

spike_threshold = -20; % mV
time_between_spikes = 5; % ms
dt = .002; %
total_time = 2000; %ms
t = 0:dt:total_time;
num_steps = length(t);
Ie_vals = linspace(106,145, 40);
firing_rates = zeros(1, 40);

%Initialization
for k = 1:40
    t0 = 40; % ms 
    I0 = Ie_vals(k); % pA
    I_e = zeros(1, num_steps);
    I_e(t >= t0) = I0;
    
    I_Na = zeros(1, num_steps);
    I_K = zeros(1, num_steps);
    I_L = zeros(1, num_steps);
    
    V = zeros(1, num_steps);
    m = zeros(1, num_steps);
    h = zeros(1, num_steps);
    n = zeros(1, num_steps);
    
    
    total_spikes = 0;
    last_spike = 0;
    all_spike_times = [];
    
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
            last_spike = t(i+1);
        end
    end
    if length(all_spike_times) > 1
        ISI_vals = diff(all_spike_times);
        if length(all_spike_times) >= 50
            firing_rates(k) = 1000/mean(ISI_vals);
        else
            firing_rates(k) = 0;
        end
    else
        firing_rates(k) = 0;
    end
end
disp(length(all_spike_times))


figure;
plot(Ie_vals, firing_rates, '-o');
xlabel('Input Current I_e (pA)');
ylabel('Firing Rate (spike/s)');
title('Current-Frequency Response (f-I Curve)');
grid on;
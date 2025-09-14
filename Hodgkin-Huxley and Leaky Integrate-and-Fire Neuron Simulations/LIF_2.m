% Model Parameters
C = 1;% nF
G_L = 50;% nS
V_L = -65;% mV
V_spk = -45;% mV
V_r = -65;% mV
t_arp = 2;% ms
I0 = 1.1; %nA

% Simulation Parameters
dt = 0.2;% ms
total_time = 2000;% ms
t = 0:dt:total_time;
num_steps = length(t);
Ie_vals = linspace(0,2, 100);
firing_rates = zeros(1, 100);


for k = 1:100
    I_e = zeros(1, num_steps);
    t0 = 40;
    I_e(t >= t0) = Ie_vals(k);
    V = zeros(1, num_steps);
    V(1) = V_L;
    refractory_time = 0;
    all_spike_times = [];
    
    for i = 1:num_steps-1
        if refractory_time > t(i)
            V(i+1) = V_r;
        else
            dV_dt = (-G_L/1000 * (V(i) - V_L) + I_e(i)) / C;
            V(i + 1) = V(i) + (dV_dt * dt);
            if V(i+1) >= V_spk
                refractory_time = t(i+1) + t_arp;
                V(i+1) = V_r;
                all_spike_times = [all_spike_times, t(i)];
            end
        end
    end
    if length(all_spike_times) > 1
        ISI_vals = diff(all_spike_times);
        firing_rates(k) = 1000/mean(ISI_vals);
    else
        firing_rates(k) = 0;
    end
end

figure;
plot(Ie_vals, firing_rates, '-o');
xlabel('Input Current I_e (nA)');
ylabel('Firing Rate (spike/s)');
title('Current-Frequency Response (f-I Curve)');
grid on;
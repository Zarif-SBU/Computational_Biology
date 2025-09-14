% Model Parameters
C = 1;% nF
G_L = 50;% nS
V_L = -65;% mV
V_spk = -45;% mV
V_r = -65;% mV
t_arp = 2;% ms
I0 = 1.1; %nA

% Simulation Parameters
dt = 0.1;% ms
total_time = 200;% ms
t = 0:dt:total_time;
num_steps = length(t);

I_e = zeros(1, num_steps);
t0 = 40;
I_e(t >= t0) = I0;
V = zeros(1, num_steps);
V(1) = V_L;
refractory_time = 0;


tic
for i = 1:num_steps-1
    if refractory_time > t(i)
        V(i+1) = V_r;
    else
        dV_dt = (-G_L * (V(i) - V_L) + I_e(i)) / (C*1000);
        V(i + 1) = V(i) + (dV_dt * dt);
        if V(i+1) >= V_spk
            refractory_time = t(i+1) + t_arp;
            V(i+1) = V_r;
        end
    end
end
toc


figure;
subplot(2,1,1)
plot(t, V);
xlabel('Time (ms)');
ylabel('V (mV)');
title('Membrane Potential');
grid on;
ylim([min(V)-5 max(V)+5]);
subplot(2,1,2)
plot(t, I_e);
xlabel('Time (ms)');
ylabel('I_e (nA)');
title('External Current');
grid on;
ylim([min(I_e) max(I_e)+.1]);
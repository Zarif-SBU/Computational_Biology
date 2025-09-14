cl_omega = 40; %cl translation rate
cl_Xp = .8; %cl protein degredation rate
cl_mu = 40; %cl transcription rate
cl_Xr = .8; %cl rna degredatation rate
cl_khalf = 15; % cl k_half
stress_factor = 0;%////////////////////////////////////// stress factor////////////////////////////////////////////////////////

cro_Khalf = 10; %cro k_half 
cro_omega = 35; %cro translation rate
cro_Xp = 1.2; % cro protein degredation rate
cro_mu = 35; % cro transcription rate
cro_Xr = 1.2; %cro rna degredation rate
time = 40; %total run time of the simulation
time_step = .01;
num_steps = time/time_step; %num of steps taken calculation from total time divided by the time step

%intializing arrays to hold all concentrations and the time
cl_pro = zeros(1, num_steps);
cl_rna = zeros(1, num_steps);
cro_pro = zeros(1, num_steps);
cro_rna = zeros(1, num_steps);
t_values = zeros(1, num_steps);

%starting concentrations
cl_rna(1) = 0;
cl_pro(1) = 0;
cro_rna(1) = 0;
cro_pro(1) = 0;

for i = 1:num_steps-1
    t_values(i+1) = t_values(i) + time_step;     %calculating t of next step

    %calculating all rate of change of concentrations
    v1 = cl_omega * cl_rna(i);
    v2 = cl_Xp * cl_pro(i);
    v2_half = stress_factor * cl_pro(i);%////////////////////////////////////// degradation rate from stress////////////////////////////////////////////////////////
    v3 = cl_mu * (1- cro_pro(i)^2/(cro_Khalf^2 + cro_pro(i)^2));
    v4 = cl_rna(i) * cl_Xr;
    v5 = cro_omega * cro_rna(i);
    v6 = cro_Xp * cro_pro(i);
    v7 = cro_mu * (1- cl_pro(i)^2 /(cl_khalf^2 + cl_pro(i)^2));
    v8 = cro_Xr * cro_rna(i);

    %Forward Eular implementation for calculating concentration at next step
    cl_pro(i+1) = cl_pro(i) + time_step * (v1 - (v2 + v2_half)); %////////////////////////////////////// concentration change from stress degradation ////////////////////////////////////////////////////////
    cl_rna(i+1) = cl_rna(i) + time_step * (v3 - v4);

    cro_pro(i+1) = cro_pro(i) + time_step * (v5 - v6);
    cro_rna(i+1) = cro_rna(i) + time_step * (v7 - v8);
end
figure;
hold on;
plot(t_values, cl_pro, 'r', 'DisplayName', "cl_{pro}, cl\_pro(1) = 0");
plot(t_values, cl_rna, 'g', 'DisplayName', "cl_{rna}, cl\_rna(1) = 0");
plot(t_values, cro_pro, 'b', 'DisplayName', "cro_{pro}, cro\_pro(1) = 0");
plot(t_values, cro_rna, 'k', 'DisplayName', "cro_{rna}, cro\_rna(1) = 0");
title('Concentration over Time');
xlabel('Time')
ylabel('Concentration')
legend;
grid on;
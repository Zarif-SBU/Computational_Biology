cl_omega = 50;
cl_Xp = 1.2;
cl_mu = 50;
cl_Xr = 1.2;
cl_khalf = 10;

cro_Khalf = 10;
cro_omega = 50;
cro_Xp = 0.8;
cro_mu = 50;
cro_Xr = 0.8;

time = 20;
time_step = .01;
num_steps = time/time_step;

figure;
title('cro_{pro} vs cl_{pro} phase plane(Eular Algorithm)');
hold on;

for j = 0:500:2000
    for k = 0:500:2000
        cl_pro = zeros(1, num_steps);
        cl_rna = zeros(1, num_steps);
        cro_pro = zeros(1, num_steps);
        cro_rna = zeros(1, num_steps);
        t_values = zeros(1, num_steps);

        cro_rna(1) = k;
        cl_rna(1) = j;
        cro_pro(1) = 0;
        cl_pro(1) = 0;

        for i = 1:num_steps-1
            t_values(i+1) = t_values(i) + time_step;
            v1 = cl_omega * cl_rna(i);
            v2 = cl_Xp * cl_pro(i);
        
            v3 = cl_mu * (1 - cro_pro(i)^2 / (cro_Khalf^2 + cro_pro(i)^2));
            v4 = cl_rna(i) * cl_Xr;
        
            v5 = cro_omega * cro_rna(i);
            v6 = cro_Xp * cro_pro(i);
        
            v7 = cro_mu * (1 - cl_pro(i)^2 / (cl_khalf^2 + cl_pro(i)^2));
            v8 = cro_Xr * cro_rna(i);
        
            cl_pro(i+1) = cl_pro(i) + time_step * (v1 - v2);
            cl_rna(i+1) = cl_rna(i) + time_step * (v3 - v4);
        
            cro_pro(i+1) = cro_pro(i) + time_step * (v5 - v6);
            cro_rna(i+1) = cro_rna(i) + time_step * (v7 - v8);
        end

        % Plot the trajectory for the current initial conditions
        plot(cl_pro, cro_pro);
    end
end

xlabel('cl_{pro}');
ylabel('cro_{pro}');
%axis([0 700 0 1300]);

hold off;
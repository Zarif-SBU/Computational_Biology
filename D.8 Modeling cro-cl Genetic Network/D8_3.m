cl_omega = 50;
cl_Xp = .8;
cl_mu = 50;
cl_Xr = .8;
cl_khalf = 15;

cro_Khalf = 10;
cro_omega = 30;
cro_Xp = 1.2;
cro_mu = 30;
cro_Xr = 1.2;

figure;
hold on;

nums_steps = 50000;
num_reactions = 8;
%for loop repeating the simulation 20 times
for k = 1:20
    %initializing concentrations
    cl_pro = zeros(1, nums_steps);
    cl_rna = zeros(1, nums_steps);
    cro_pro = zeros(1, nums_steps);
    cro_rna = zeros(1, nums_steps);
    t_values = zeros(1, nums_steps);
    cl_rna(1) = 0;
    cl_pro(1) = 0;
    cro_rna(1) = 0;
    cro_pro(1) = 0;
    for i= 1:nums_steps-1
        %calculating each rection rate
        rxn_rate = zeros(1, num_reactions);
        rxn_rate(1) = cl_omega * cl_rna(i);
        rxn_rate(2) = cl_Xp * cl_pro(i);
    
        rxn_rate(3) = cl_mu * (1- cro_pro(i)^2/(cro_Khalf^2 + cro_pro(i)^2));
        rxn_rate(4) = cl_rna(i) * cl_Xr;
    
        rxn_rate(5) = cro_omega * cro_rna(i);
        rxn_rate(6) = cro_Xp * cro_pro(i);
    
        rxn_rate(7) = cro_mu * (1- cl_pro(i)^2 /(cl_khalf^2 + cl_pro(i)^2));
        rxn_rate(8) = cro_Xr * cro_rna(i);
        %calculate next time step
        alpha0 = sum(rxn_rate);
        y = rand();
        tau = -log(y)/alpha0;
        t_values(i+1) = t_values(i) + tau;
        %calculate next reaction
        y = rand();
        y_mod = alpha0*y;
        next_rxn = 0;
        for j = 1:num_reactions
            if y_mod < sum(rxn_rate(1:j))
                next_rxn = j;
                break;
            end
        end
        %select next reaction
        if (next_rxn == 1)
            cl_pro(i+1) = cl_pro(i) + 1;
            cl_rna(i+1) = cl_rna(i);
            cro_pro(i+1) = cro_pro(i);
            cro_rna(i+1) = cro_rna(i);
        elseif (next_rxn == 2)
            cl_pro(i+1) = cl_pro(i) - 1;
            cl_rna(i+1) = cl_rna(i);
            cro_pro(i+1) = cro_pro(i);
            cro_rna(i+1) = cro_rna(i);
        elseif (next_rxn == 3)
            cl_pro(i+1) = cl_pro(i);
            cl_rna(i+1) = cl_rna(i) + 1;
            cro_pro(i+1) = cro_pro(i);
            cro_rna(i+1) = cro_rna(i);
        elseif (next_rxn == 4)
            cl_pro(i+1) = cl_pro(i);
            cl_rna(i+1) = cl_rna(i) - 1;
            cro_pro(i+1) = cro_pro(i);
            cro_rna(i+1) = cro_rna(i);
        elseif (next_rxn == 5)
            cl_pro(i+1) = cl_pro(i);
            cl_rna(i+1) = cl_rna(i);
            cro_pro(i+1) = cro_pro(i) + 1;
            cro_rna(i+1) = cro_rna(i);
        elseif (next_rxn == 6)
            cl_pro(i+1) = cl_pro(i);
            cl_rna(i+1) = cl_rna(i);
            cro_pro(i+1) = cro_pro(i) - 1;
            cro_rna(i+1) = cro_rna(i);
        elseif (next_rxn == 7)
            cl_pro(i+1) = cl_pro(i);
            cl_rna(i+1) = cl_rna(i);
            cro_pro(i+1) = cro_pro(i);
            cro_rna(i+1) = cro_rna(i) + 1;
        elseif (next_rxn == 8)
            cl_pro(i+1) = cl_pro(i);
            cl_rna(i+1) = cl_rna(i);
            cro_pro(i+1) = cro_pro(i);
            cro_rna(i+1) = cro_rna(i) - 1;
        end
    end
    subplot(2,1,1); % Concentration vs Time
    hold on;
    stairs(t_values, cl_pro, 'r');
    stairs(t_values, cl_rna, 'g');
    stairs(t_values, cro_pro, 'b');
    stairs(t_values, cro_rna, 'k');
    legend('cl_{pro}', 'cl_{rna}', 'cro_{pro}', 'cro_{rna}')
    xlabel('Time');
    ylabel('Concentration');
    grid on;
    title('All Concentrations');

    subplot(2,1,2); % Cro_pro vs Cl_pro
    stairs(cl_pro, cro_pro);
    hold on;
    xlabel('Cl_{Pro}');
    ylabel('Cro_{Pro}');
    grid on;
    title('cl_{pro} vs cro_{pro}');
    axis([0 600 0 1000])
end
%legend('cl_{pro}', 'cl_{rna}', 'cro_{pro}', 'cro_{rna}')
%xlabel('Time');
%ylabel('Concentration');
%grid on;
hold off;
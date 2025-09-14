K_half = .33; %mM
mu = 1; %s-1
Xr = 1; %s-1
Xp = 1; %s-1 
omega = 1; %s-1
Dr = 1e-4; %um * s^-1
Dp = 1e-4; %um * s^-1

%Setting time and space constraints
total_time = 120; %s
time_step = .01; %s
num_time_steps = total_time/time_step;
total_dist = 3.0; %um
dist_step = .02; %um
num_dist_steps = ceil(total_dist/dist_step)+1;

% Setting our arrays to store discritized time and space values
t_values = 0:time_step:total_time;
r_values = 0:dist_step:total_dist;

% 2d arrays storing the concentration of rna and protein respective to time
% and space
pro = zeros(num_dist_steps, num_time_steps+1);
rna = zeros(num_dist_steps, num_time_steps+1);

%Initial conditions
pro(1, 1) = 1;
rna(1, 1) = 1;

% Calculating rate of change and applying to the Forward Eular method
for i = 1:num_time_steps % Going over time
    for j = 1:num_dist_steps % Going over our 1d space
        if j == 1 %Calculating our laplacian operator
            laplacian_rna = (rna(2, i) - 2 *rna(1,i))/ (dist_step^2);
            laplacian_pro = (pro(2, i) - 2 *pro(1,i))/ (dist_step^2); 
        elseif j == num_dist_steps
            laplacian_rna = (rna(j-1, i) - 2 *rna(j,i))/ (dist_step^2);
            laplacian_pro = (pro(j-1, i) - 2 *pro(j,i))/ (dist_step^2); 
        else
            laplacian_rna = (rna(j+1, i) + rna(j-1, i) -2 * rna(j, i))/(dist_step^2);
            laplacian_pro = (pro(j+1, i) + pro(j-1, i) -2 * pro(j, i))/(dist_step^2);
        end
        % Rate of change respective to both time and change
        dr_dt = (mu * pro(j,i)^2)/(K_half^2 + pro(j,i)^2) - Xr * rna(j,i) + Dr * laplacian_rna;
        dp_dt = omega * rna(j,i) - Xp * pro(j,i) + Dp * laplacian_pro;
        % Forward Eular Method
        rna(j, i+1) = rna(j,i) + time_step * dr_dt;
        pro(j, i+1) = pro(j,i) + time_step * dp_dt;
    end
end


% Plotting results
figure;
subplot(3,1,1);
hold on;
plot(t_values, rna(1, :));
plot(t_values, pro(1, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration at First Position');
legend('RNA', 'Protein');
grid on;

subplot(3,1,2);
hold on;
plot(t_values, rna(2, :));
plot(t_values, pro(2, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration at Second Position');
legend('RNA', 'Protein');
grid on;

subplot(3,1,3);
hold on;
plot(t_values, rna(3, :));
plot(t_values, pro(3, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration at Third Position');
legend('RNA', 'Protein');
grid on;

figure;

subplot(2,1,1);
hold on;
plot(r_values, rna(:, 1), 'b');
plot(r_values, pro(:, 1), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Initial RNA and Protein Concentration vs Distance');
legend('RNA', 'Protein');
grid on;

subplot(2,1,2);
hold on;
plot(r_values, rna(:, end), 'b');
plot(r_values, pro(:, end), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Final RNA and Protein Concentration vs Distance');
legend('RNA', 'Protein');
grid on;
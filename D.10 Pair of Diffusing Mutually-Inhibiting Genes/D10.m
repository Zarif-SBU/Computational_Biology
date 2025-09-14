% Setting model paramenters
K_half = 0.33; % mM
mu = 1; % s^-1
chi_r = 1; % s^-1
chi_p = 1; % s^-1 
omega = 1; % s^-1
Dr = 1e-4; % um^2/s
Dp = 1e-4; % um^2/s

% Setting simulation parameters
total_time = 6; % s
time_step = 0.01; % s
num_time_steps = total_time / time_step;
total_dist = 3.0; % um
dist_step = 0.02; % um
num_dist_steps = total_dist / dist_step + 1;
t_values = 0:time_step:total_time;
r_values = 0:dist_step:total_dist;

% Initialize concentration matrices
X_rna = zeros(num_dist_steps, num_time_steps+1);
X_pro = zeros(num_dist_steps, num_time_steps+1);
Y_rna = zeros(num_dist_steps, num_time_steps+1);
Y_pro = zeros(num_dist_steps, num_time_steps+1);
% Initial conditions
X_rna(1, 1) = 1;
X_pro(1, 1) = 1;
Y_rna(end, 1) = 1;
Y_pro(end, 1) = 1;

% Calculating rate of change and applying to the Forward Eular method
for i = 1:num_time_steps % Going over time
    for j = 1:num_dist_steps % Going over our 1d space
        if j == 1 %Calculating our laplacian operator
            laplacian_X_rna = (X_rna(2, i) - 2 *X_rna(1,i))/ (dist_step^2);
            laplacian_X_pro = (X_pro(2, i) - 2 *X_pro(1,i))/ (dist_step^2); 
            laplacian_Y_rna = (Y_rna(2, i) - 2 *Y_rna(1,i))/ (dist_step^2);
            laplacian_Y_pro = (Y_pro(2, i) - 2 *Y_pro(1,i))/ (dist_step^2); 
        elseif j == num_dist_steps
            laplacian_X_rna = (X_rna(j-1, i) - 2 *X_rna(j,i))/ (dist_step^2);
            laplacian_X_pro = (X_pro(j-1, i) - 2 *X_pro(j,i))/ (dist_step^2); 
            laplacian_Y_rna = (Y_rna(j-1, i) - 2 *Y_rna(j,i))/ (dist_step^2);
            laplacian_Y_pro = (Y_pro(j-1, i) - 2 *Y_pro(j,i))/ (dist_step^2); 
        else
            laplacian_X_rna = (X_rna(j+1, i) + X_rna(j-1, i) -2 * X_rna(j, i))/(dist_step^2);
            laplacian_X_pro = (X_pro(j+1, i) + X_pro(j-1, i) -2 * X_pro(j, i))/(dist_step^2);
            laplacian_Y_rna = (Y_rna(j+1, i) + Y_rna(j-1, i) -2 * Y_rna(j, i))/(dist_step^2);
            laplacian_Y_pro = (Y_pro(j+1, i) + Y_pro(j-1, i) -2 * Y_pro(j, i))/(dist_step^2);
        end

        % Calculating rate of change respective to both time and change
        dXrna_dt = mu * (1 - (Y_pro(j, i)^2 / (K_half^2 + Y_pro(j, i)^2))) - chi_r * X_rna(j, i) + Dr * laplacian_X_rna;
        dXpro_dt = omega * X_rna(j, i) - chi_p * X_pro(j, i) + Dp * laplacian_X_pro;
        dYrna_dt = mu * (1 - (X_pro(j, i)^2 / (K_half^2 + X_pro(j, i)^2))) - chi_r * Y_rna(j, i) + Dr * laplacian_Y_rna;
        dYpro_dt = omega * Y_rna(j, i) - chi_p * Y_pro(j, i) + Dp * laplacian_Y_pro;

        % Forward Eular Method
        X_rna(j, i+1) = X_rna(j, i) + time_step * dXrna_dt;
        X_pro(j, i+1) = X_pro(j, i) + time_step * dXpro_dt;
        Y_rna(j, i+1) = Y_rna(j, i) + time_step * dYrna_dt;
        Y_pro(j, i+1) = Y_pro(j, i) + time_step * dYpro_dt;
    end
end

% Plotting results
figure;
subplot(3,1,1);
hold on;
plot(t_values, X_rna(1, :));
plot(t_values, X_pro(1, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration of Species X at First Position');
legend('RNA', 'Protein');
grid on;

subplot(3,1,2);
hold on;
plot(t_values, X_rna(ceil(num_dist_steps/2), :));
plot(t_values, X_pro(ceil(num_dist_steps/2), :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration of Species X at Middle Position');
legend('RNA', 'Protein');
grid on;

subplot(3,1,3);
hold on;
plot(t_values, X_rna(end, :));
plot(t_values, X_pro(end, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration of Species X at Last Position');
legend('RNA', 'Protein');
grid on;

figure;
subplot(3,1,1);
hold on;
plot(t_values, Y_rna(1, :));
plot(t_values, Y_pro(1, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration of Species Y at First Position');
legend('RNA', 'Protein');
grid on;

subplot(3,1,2);
hold on;
plot(t_values, Y_rna(ceil(num_dist_steps/2), :));
plot(t_values, Y_pro(ceil(num_dist_steps/2), :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration of Species Y at Middle Position');
legend('RNA', 'Protein');
grid on;

subplot(3,1,3);
hold on;
plot(t_values, Y_rna(end, :));
plot(t_values, Y_pro(end, :));
xlabel('Time (s)');
ylabel('Concentration (mM)');
title('RNA and Protein Concentration of Species Y at Last Position');
legend('RNA', 'Protein');
grid on;

figure;
subplot(3,1,1);
hold on;
plot(r_values, X_rna(:, 1), 'b');
plot(r_values, X_pro(:, 1), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Initial Concentration of Species X at all distances');
legend('RNA', 'Protein');
grid on;

subplot(3,1,2);
hold on;
plot(r_values, X_rna(:, ceil(num_dist_steps/2)), 'b');
plot(r_values, X_pro(:, ceil(num_dist_steps/2)), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Middle Concentration of Species X at all distances');
legend('RNA', 'Protein');
grid on;

subplot(3,1,3);
hold on;
plot(r_values, X_rna(:, end), 'b');
plot(r_values, X_pro(:, end), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Final Concentration of Species X at all distances');
legend('RNA', 'Protein');
grid on;

figure;
subplot(3,1,1);
hold on;
plot(r_values, Y_rna(:, 1), 'b');
plot(r_values, Y_pro(:, 1), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Initial Concentration of Species Y at all distances');
legend('RNA', 'Protein');
grid on;

subplot(3,1,2);
hold on;
plot(r_values, Y_rna(:, ceil(num_dist_steps/2)), 'b');
plot(r_values, Y_pro(:, ceil(num_dist_steps/2)), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Middle Concentration of Species Y at all distances');
legend('RNA', 'Protein');
grid on;

subplot(3,1,3);
hold on;
plot(r_values, Y_rna(:, end), 'b');
plot(r_values, Y_pro(:, end), 'r');
xlabel('Distance (µm)');
ylabel('Concentration (mM)');
title('Final Concentration of Species Y at all distances');
legend('RNA', 'Protein');
grid on;
%D7
%{
omega = 1; % Protein synthesis rate 
Xr = 1; % Rna degredation rate
Xp = 1; % Protein degredation rate
mu = 2; %rna synthesis
h = 2;
k_half = 1;

init_pro = 0;
init_rna = 0;

step_size = 0.01; % Step size
time = 20; % total time in seconds
num_steps = time/step_size;
Pro = zeros(1, num_steps); % Protein amounts
Rna = zeros(1, num_steps); % RNA amounts
t_values = zeros(1, num_steps);

Pro(1) = init_pro;
Rna(1) = init_rna;
t = 0;

for i = 1:num_steps-1
    t = t + step_size;
    dP = omega * Rna(i) - Xp * Pro(i);
    dR = (mu * Pro(i) ^ h /( k_half ^ h + Pro(i) ^ h)) - Xr * Rna(i);
    Pro(i+1) = Pro(i) + dP * step_size;
    Rna(i+1) = Rna(i) + dR * step_size;
    t_values(i+1) = t; %we i + 1 to skip first index since we know that will 0
end

figure;
hold on;
plot(t_values, Pro, 'DisplayName', 'Protein');
plot(t_values, Rna, 'DisplayName', 'RNA');
xlabel('Time (s)');
ylabel('Concentration (mM)');
legend;
title('Protein and RNA Dynamics');
hold off;
%}
%{
omega = 1; % Protein synthesis rate 
Xr = 1; % Rna degredation rate
Xp = 1; % Protein degredation rate
mu = 1; %rna synthesis
h = 2;
k_half = .33;

step_size = 0.001; % Step size
time = 20; % total time in seconds
num_steps = time/step_size;

        figure;
        hold on;
colors = jet(64); % Generate a colormap with 64 colors
counter = 1;
for j = 0: .2:1.4
    for k = 0: .2: 1.4
        Pro = zeros(1, num_steps + 1); % Protein amounts
        Rna = zeros(1, num_steps + 1); % RNA amounts
        Pro(1) = j;
        Rna(1) = k;
        for i = 1:num_steps
                dP = omega * Rna(i) - Xp * Pro(i);
                dR = (mu * Pro(i) ^ h /( k_half ^ h + Pro(i) ^ h)) - Xr * Rna(i);
                Pro(i+1) = Pro(i) + dP * step_size;
                Rna(i+1) = Rna(i) + dR * step_size;
        end
        plot(Pro, Rna, 'Color', colors(counter, :));
        plot(Pro(num_steps+1), Rna(num_steps+1), 'o', 'LineWidth', 3, 'Color', 'blue');
        counter = counter + 1;
    end
end
xlabel('Protein concentration (mM)');
ylabel('Rna concentraion (mM)');
title('Protein and RNA Dynamics');
hold off;

%}

omega = 1; % Protein synthesis rate 
Xr = 1; % Rna degredation rate
Xp = 1; % Protein degredation rate
mu = 2; %rna synthesis
h = 2;
k_half = 1;

Pro = linspace(0, 2, 100);

Rna_nullcline = (mu * Pro.^h) ./ (Xr * (k_half^h + Pro.^h));

Pro_nullcline = (Xp / omega) * Pro;

syms Xp_sym
eq = (mu * Xp_sym^h) / (Xr * (k_half^h + Xp_sym^h)) == (Xp / omega) * Xp_sym;
Xp_eq = double(vpasolve(eq, Xp_sym)); 
Xr_eq = (Xp / omega) * Xp_eq;

figure;
hold on;
plot(Pro, Rna_nullcline, 'DisplayName', 'RNA Nullcline');
plot(Pro, Pro_nullcline, 'DisplayName', 'Protein Nullcline');
plot(Xp_eq, Xr_eq, 'o', 'DisplayName', 'Equilibrium Points');

xlabel('Protein concentration [X_p]');
ylabel('RNA concentration [X_r]');
title('Nullclines and Stationary Points');
legend;
hold off;
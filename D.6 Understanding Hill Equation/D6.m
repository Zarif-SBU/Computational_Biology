%D6
%Vary h
%{
V_max = 5;
K_half = 20;
h = 1;
h_two = 2;
h_three = 10;

S_values = linspace(0,100,101);
v_values = zeros(1, 101);
v_values_two = zeros(1, 101);
v_values_three = zeros(1, 101);

function v = hill_equation(V_max, h, K_half, S)
    v = V_max * S ^ h / (K_half + S ^ h);
end

for i = 1:101
    v_values(i) = hill_equation(V_max, h, K_half, S_values(i));
    v_values_two(i) = hill_equation(V_max, h_two, K_half, S_values(i));
    v_values_three(i) = hill_equation(V_max, h_three, K_half, S_values(i));
end

figure;
hold on;
plot(S_values, v_values, 'b');
plot(S_values, v_values_two, 'r');
plot(S_values, v_values_three, 'g');
xlabel('Substrate concentration(mM)');
ylabel('Rate(mM/s)');
title('Vary h');
legend('h = 1', 'h = 2', 'h = 10');
%}



%Vary K_half
%{
V_max = 5;
K_half = 10;
K_half_two = 20;
K_half_three = 40;
h = 2;

S_values = linspace(0,100,101);
disp(S_values)
v_values = zeros(1, 101);
v_values_two = zeros(1, 101);
v_values_three = zeros(1, 101);

function v = hill_equation(V_max, h, K_half, S)
    v = V_max * S ^ h / (K_half + S ^ h)
end

for i = 1:101
    v_values(i) = hill_equation(V_max, h, K_half, S_values(i));
    v_values_two(i) = hill_equation(V_max, h, K_half_two, S_values(i));
    v_values_three(i) = hill_equation(V_max, h, K_half_three, S_values(i));
end

figure;
hold on;
plot(S_values, v_values, 'b');
plot(S_values, v_values_two, 'r');
plot(S_values, v_values_three, 'g');
xlabel('Substrate concentration(mM)');
ylabel('Rate(mM/s)');
title('Vary K_{half}');
legend('k_{1/2} = 10mM' , 'k_{1/2} = 20mM', 'k_{1/2} = 40mM');
%}


%Vary V_max
%{
V_max = 2;
V_max_two = 5;
V_max_three = 10;
K_half = 20;
h = 2;

S_values = linspace(0,100,101);
v_values = zeros(1, 101);
v_values_two = zeros(1, 101);
v_values_three = zeros(1, 101);

function v = hill_equation(V_max, h, K_half, S)
    v = V_max * S ^ h / (K_half + S ^ h);
end

for i = 1:101
    v_values(i) = hill_equation(V_max, h, K_half, S_values(i));
    v_values_two(i) = hill_equation(V_max_two, h, K_half, S_values(i));
    v_values_three(i) = hill_equation(V_max_three, h, K_half, S_values(i));
end

figure;
hold on;
plot(S_values, v_values, 'b');
plot(S_values, v_values_two, 'r');
plot(S_values, v_values_three, 'g');
xlabel('Substrate concentration(mM)');
ylabel('Rate(mM/s)');
title('Vary V_{max}');
legend('V_{max} = 2mM/s', 'V_{max} = 5mM/s', 'V_{max} = 10mM/s');
%}

%for a Hill equation with h = 4, V_max = 20, K_half = 50, the graph would
%have a ceiling at 20 and it will reach the celing reletively quickly,
%maybe by the time it reaches a substrate concentration of 22.5
%{
V_max = 20;
K_half = 50;
h = 4;

S_values = linspace(0,100,101);
v_values = zeros(1, 101);
v_values_two = zeros(1, 101);
v_values_three = zeros(1, 101);

function v = hill_equation(V_max, h, K_half, S)
    v = V_max * S ^ h / (K_half + S ^ h);
end

for i = 1:101
    v_values(i) = hill_equation(V_max, h, K_half, S_values(i));
end

figure;
hold on;
plot(S_values, v_values, 'b');
xlabel('Substrate concentration(mM)');
ylabel('Rate(mM/s)');
title('Rate  at substrate concentration');
legend('V_{max} = 20, h = 2, K_{half} = 50');
%}

%question 3
V_max = 20;
K_half = 50;
h = 4;


init_S = 100;
step_size = 0.01; % Step size
num_steps = 10/step_size;
S_values = zeros(1, num_steps + 1); % Protein amounts
S_values(1) = init_S;
t_values = zeros(1, num_steps + 1);
t = 0;

function v = hill_equation(V_max, h, K_half, S)
    v = V_max * S ^ h / (K_half + S ^ h);
end


for i = 1:num_steps
    t = t + step_size;
    S_values(i+1) = S_values(i) + hill_equation(V_max, h, K_half, S_values(i)) * step_size;
    t_values(i+1) = t;
end

figure;
hold on;
plot(t_values, S_values, 'b');
xlabel('Time(sec)');
ylabel('Substrate concentration(mM)');
title('Substrate concentration over time');
legend('V_{max} = 20, h = 2, K_{half} = 50');



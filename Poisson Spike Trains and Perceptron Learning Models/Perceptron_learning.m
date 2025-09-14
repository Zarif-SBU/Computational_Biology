n = 1;
N = 2;
w = zeros(1,N+1);
T = 1000; % Stop condition

patterns = [-1 -1 -1; 1 -1 -1; -1 1 -1; 1 1 -1];

AND = [-1 -1 -1 1];
OR = [-1 1 1 1];
XOR = [-1 1 1 -1];
consective_correct = 0;
convergence_time = T;

performance = zeros(1, T);
for i = 1:T
    r = randi(4);
    x = patterns(r, :);
    y_target = XOR(r); %% Change this to AND OR or XOR depending on what you are looking for
    y_actual = sign(w*x');
    if y_actual == y_target
        performance(i) = 1;
        consective_correct = consective_correct + 1;
        if consective_correct == 200
            convergence_time = i-200;
        end
    else
        consective_correct = 0;
    end
    for j = 1:N+1
        w(j) = w(j) + n * (y_target - y_actual) * x(j);
    end
end

disp(convergence_time)
figure('Position', [600 400 500 200])
plot(performance, '.k', 'MarkerSize', 1)
xlabel('Step#')
ylabel('Performance(fail: 0, pass: 1)')
title('Performance each step for XOR')
ylim([-0.2 1.2])

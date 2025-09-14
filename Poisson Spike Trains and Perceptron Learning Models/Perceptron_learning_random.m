n = 1;
N = 50;
M = 40; %edit number of patters
T = 1000; % Stop condition
runs = 1; %edit number of runs
convergences = zeros(1, runs);

figure('Position', [600 400 500 1000])

for k = 1:runs
    w = zeros(1,N+1);
    patterns = randi(2,N,M)-1;
    patterns = 2 * patterns - 1;
    yt = [ones(1, M/2), -ones(1, M/2)];
    yt = yt(randperm(M));
    consective_correct = 0;

    convergence_time = T;
    performance = zeros(1, T);

    for i = 1:T
        r = randi(M);
        x = [patterns(:, r)' 1];
        y_target = yt(r); %% Change this to AND OR or XOR depending on what you are looking for
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
    
    convergences(1, k) = convergence_time;
    if k <=5
        subplot(5, 1, k)
        plot(performance, '.k', 'MarkerSize', 1)
        xlabel('Step#')
        ylabel('Performance (0 or 1)')
        title('Performance each step')
        ylim([-0.2 1.2])
    end
end
converged_runs = convergences < T;
average_convergence_time = mean(convergences(converged_runs))/M;
num_not_converged = sum(~converged_runs);
fprintf('Average convergence time per pattern (successful runs): %.2f steps\n', average_convergence_time);
fprintf('Number of runs that did NOT converge: %d\n', num_not_converged);
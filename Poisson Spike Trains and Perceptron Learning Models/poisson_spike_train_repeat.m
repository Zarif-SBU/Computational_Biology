N = 50; % Spike train count
lambda = 10; % spikes/s
T_vals = [10,100,200,400,800,1600,3200,6400]; % Length of each

CVs = zeros(1, length(T_vals));
FFs = zeros(1, length(T_vals));

%% Get the spike trains for 50 trials for each time value
for t_i = 1:length(T_vals)
    T = T_vals(t_i);
    max_spikes = 2*lambda*T;
    spike_trains = NaN(N, max_spikes);

    %% Calculate spike trains for current simulation time
    for i = 1:N
        spike_count = 0;
        current_time = 0;
        while true
            u = rand(1);
            ISI = -log(u)/lambda;
            current_time = current_time + ISI;
            if current_time > T
                break
            end
            spike_count = spike_count + 1;
            spike_trains(i, spike_count) = current_time;
        end
    end

    %% Calculate CV of Trial 1
    spike_count_trial_one = spike_trains(1, ~isnan(spike_trains(1, :)));
    ISI = diff(spike_count_trial_one);
    CVs(t_i) = std(ISI)/mean(ISI);

    %% Calculate number of spikes in each trial to calculate variance between trial for FF
    spike_counts = sum(~isnan(spike_trains), 2)';
    FFs(t_i) = var(spike_counts)/mean(spike_counts);
end

figure;
semilogx(T_vals, CVs, '-o','DisplayName', 'CV')
hold on;
semilogx(T_vals, FFs, '-s', 'DisplayName', 'FF')
xlabel('Time (sec)')
ylabel('Value')
title('CV and FF at different simulation length')
grid on;
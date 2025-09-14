N = 50; % Spike train count
lambda = 10; % spikes/s
T = 5; % Length of each
max_spikes = 80;

spike_trains = NaN(N, max_spikes);

u = rand(100,1);

%% Get the spike trains for 50 trials for which we will create a rastor plo
for i = 1:N
    spike_count = 0;
    current_time = 0;
    while true
        u = rand(1);
        ISI = -log(u)/lambda;
        current_time = current_time + ISI;
        if spike_count >= max_spikes || current_time > T
            break
        end
    spike_count = spike_count + 1;
    spike_trains(i, spike_count) = current_time;
    end 
end

%% Get Average Spikes
avg_spikes = sum(~isnan(spike_trains), "all")/(N * T);
disp(avg_spikes)

%% build the PSTH data
bin_count = T/.2;
PSTH = zeros(1, bin_count);
bin_edges = 0:.2:5;
bin_middles = bin_edges(1:end-1) + .1;

for i = 1:bin_count
    current_time = bin_edges(i);
    next_time = bin_edges(i + 1);
    % Count spikes in this bin across all trials
    total_spikes = 0;
    for trial_num = 1:N
        % Get spikes for this trial (excluding NaNs)
        spikes = spike_trains(trial_num, ~isnan(spike_trains(trial_num, :)));
        % Count spikes in [current_time, next_time)
        spikes_in_bin = (spikes >= current_time) & (spikes < next_time);
        total_spikes = total_spikes + sum(spikes_in_bin);
    end

    PSTH(i) = total_spikes / (N * .2);
end

%% Calculate the CV of the ISIs in each trial
spike_counts = zeros(1, N);
CVs = zeros(1, N);
for i = 1:N
    spikes = spike_trains(i, ~isnan(spike_trains(i, :)));
    spike_counts(i) = length(spikes);
    ISI = diff(spikes);
    mean_ISI = mean(ISI);
    std_ISI = std(ISI);
    CVs(i) = std_ISI/mean_ISI;
  
end

%% Calculate average CV
CV_mean = sum(CVs)/N;
disp(CV_mean)

%% Calculate Fano Factor(FF)
FF = var(spike_counts)/mean(spike_counts);
disp(FF)

figure('Position', [500 100 600 700]);
subplot(2,1,1)
plot(spike_trains, 1:N, '.k');
xlabel('Time (sec)')
ylabel('Trial #')
title('Raster Plot')
grid on;


subplot(2,1,2)
bar(bin_middles, PSTH, 'FaceColor', 'none', LineWidth=1.5);
xlabel('Time (sec)')
ylabel('Firing rate (spk/s)')
title('PSTH')
grid on;

figure;
bar(1:N, CVs)
xlabel('Trial #')
ylabel('Coefficient of variability')
title('CV plot')
grid on;
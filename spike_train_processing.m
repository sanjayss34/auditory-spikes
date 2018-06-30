function [] = spike_train_processing(trial_length, bin_length)
% This function takes CohenNeurons.mat and pre-processes the data into an
% m-by-k-by-n matrix with m matrices (1 for each neuron), k rows of trials 
% and n = 2001 columns for (1) TNR and (2000) times in spike train.
tic;
% Default bin length is 1 ms
if ~exist('bin_length', 'var')
    bin_length = 1;
end

% Default trial length is 2000ms
if ~exist('trial_length', 'var')
    trial_length = 2000;
end

% Load CohenNeurons.mat
CohenNeurons = load('CohenNeurons.mat');
CohenNeurons = CohenNeurons.CohenNeurons;

% Number of neurons
num_neurons = 16;
num_trials = 464;

% Pre-define binranges from 0 to 2000ms
binranges = linspace(0,trial_length,trial_length/bin_length+1);

% Preallocate mxkxn matrix with m neurons, k trials and n = 2000ms
% First column is TNR
% Next 2001 columns is spike train (we will remove 2001st column at the end)
spike_array = zeros(num_neurons, num_trials, trial_length/bin_length+2);

% For each neuron...
for i = 1:num_neurons
    % For each trial...
    for j = 1:num_trials
        % Create a spike train and store in spike_array
        spike_array(i,j,2:end) = histc(CohenNeurons(i).trials(j).spikes, binranges);
        % Record TNR
        spike_array(i,j,1) = CohenNeurons(i).trials(j).TNR;
    end
    % Print i
%     i
end
toc

% amplitude = stimulus_binning(trial_length/bin_length);

filtered_stimulus_initial = filter_stimulus(binranges);
filtered_stimulus = zeros(size(spike_array, 1), size(spike_array, 2), size(spike_array, 3)-1);
for i=1:size(filtered_stimulus, 1)
    for j=1:size(filtered_stimulus, 2)
        for k=1:size(filtered_stimulus, 3)
            filtered_stimulus(i, j, k) = filtered_stimulus_initial(i, j, ceil(size(filtered_stimulus_initial, 3)*k/size(filtered_stimulus, 3)));
        end
    end
end

% Remove last column of spike train (unwanted bin)
spike_array(:,:,end) = [];
filtered_stimulus(:,:,end) =[];

% Save data
save('spike_trains.mat', 'spike_array');
save('filtered_stimulus.mat', 'filtered_stimulus');
% save('amplitude.mat', 'amplitude');

% Check if any bin has more than 1 spike
s = spike_array(:,:,2:end);
% max(max(max(s)))

% figure;
% plot(reshape(filtered_stimulus(1,1,:)*10^4, [1 100]))
% hold on
% spike_array = spike_array(:,:,2:end);
% plot(reshape(spike_array(1,1,:), [1 100]))

end
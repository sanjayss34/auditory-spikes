% This file indiscriminately concatenates all trials for a single neuron.
% The result is a 16-cell array with only 1 concatenated train for each
% neuron. The trials for each neuron are sorted by increasing TNR
% intensity. In particular, for any neuron, there is 1 spike train, where
% the earliest times are from trials with TNR intensity 0 and the latest
% times are from trials with TNR intensity 85 (which is the highest intensity).

% Load sorted trains
sorted_trains = load('sorted_trains.mat');
sorted_trains = sorted_trains.sorted_trains;
load('sorted_stimuli.mat');

num_neurons = 16;
num_tnrs = 7;

% Preallocate a 16-cell array with 1 train for each neuron
neuron_trains = cell(16,1);
feature_intensity = cell(16,1);

% For each neuron...
for i = 1:num_neurons
    % Create empty vector for concatenated train
    A = [];
    % For each TNR intensity...
    for j = 1:num_tnrs
        % Concatenate the train associated with that intensity
        A = [A sorted_trains{i,j}];
    end
    neuron_trains{i,1} = A;
    B = [];
    for j=1:num_tnrs
        B = [B sorted_stimuli{i,j}];
    end
    feature_intensity{i,1} = B;
    % Print i
    i
end

save('neuron_trains.mat', 'neuron_trains');
save('feature_intensity.mat', 'feature_intensity');
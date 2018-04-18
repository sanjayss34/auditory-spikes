% This file sorts the spike trains by TNR intensity and then concatenates
% them accordingly. The output is a 16x7 cell array with 16 rows of neurons
% and 7 columns of TNR intensities. Inside each cell is the concatenated
% spike train corresponding to all trials for a specific neuron with a
% specific TNR intensity.

% Load spike trains
spike_trains = load('spike_trains.mat');
load('filtered_stimulus.mat');

spike_trains = spike_trains.spike_array;
% Remove first column (TNR intensity)
spike_trains(:,:,1) = [];

% All bins with spikes become +1
spike_trains(spike_trains > 0) = 1;
% All bins with spikes become -1
spike_trains(spike_trains == 0) = -1;

% Load TNR indices
TNR_index = load('TNR_index.mat');
TNR_index = TNR_index.ind;

% Number of TNR intensities, from 0 to 85
num_TNRs = 7;
num_neurons = 16;

% Preallocate cell array with 7 columns of TNR intensities
TNRsorted_trains = cell(1,7);
TNRstimuli = cell(1,7);

% For each TNR intensity...
for i = 1:num_TNRs
    % Find trials with TNR intensity i
    ind = (TNR_index == i);
    % Store trials with TNR intensity i in cell column i
    TNRsorted_trains{1,i} = spike_trains(:,ind,:);
    TNRstimuli{1,i} = filtered_stimulus(:,ind,:);
end

% Preallocate cell array with 16 rows of neurons and 7 columns of TNR intensities
sorted_trains = cell(16,7);
sorted_stimuli = cell(16,7);

% For each TNR intensity...
for i = 1:num_TNRs
    % Get TNR trials for that intensity
    TNRtrials = TNRsorted_trains{1,i};
    TNRtrials_stimuli = TNRstimuli{1,i};
    % For each neuron...
    for j = 1:num_neurons
        % Get specific trials for that neuron
        neurontrials = TNRtrials(j,:,:);
        s2 = size(neurontrials,2);
        s3 = size(neurontrials,3);
        neurontrials = reshape(neurontrials, [s2, s3]);
        % Concatenate trains into 1 row vector
        neurontrials = reshape(neurontrials, 1, []);
        % Save concatenated train in cell array
        sorted_trains{j,i} = neurontrials;
        % Get specific trials for that neuron
        stimulustrials = TNRtrials_stimuli(j,:,:);
        s2 = size(stimulustrials,2);
        s3 = size(stimulustrials,3);
        stimulustrials = reshape(stimulustrials, [s2, s3]);
        % Concatenate trains into 1 row vector
        stimulustrials = reshape(stimulustrials, 1, []);
        % Save concatenated train in cell array
        sorted_stimuli{j,i} = stimulustrials;
    end
    % Print i
    i
end

save('sorted_trains.mat', 'sorted_trains');
save('sorted_stimuli.mat', 'sorted_stimuli');
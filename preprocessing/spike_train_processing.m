% This function takes CohenNeurons.mat and pre-processes the data into an
% m-by-k-by-n matrix with m matrices (1 for each neuron), k rows of trials 
% and n = 2001 columns for (1) TNR and (2000) times in spike train.

% Load CohenNeurons.mat
CohenNeurons = load('CohenNeurons.mat');
CohenNeurons = CohenNeurons.CohenNeurons;

% Number of neurons
num_neurons = 16;
num_trials = 464;

% Preallocate mxkxn matrix with m neurons, k trials and n = 2000ms
% First column is TNR
% Next 2001 columns is spike train (we will remove 2001st column at the end)
spike_array = zeros(num_neurons, num_trials, 2002);

% Pre-define binranges from 0 to 2000ms
binranges = linspace(0,2000,2001);

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
    i
end

% Remove 2001st column of spike train (2002nd column of matrix)
spike_array(:,:,end) = [];

% Save data
save('spike_trains.mat', 'spike_array');
% This file produces 2 outputs. The first uses neuron_trains.mat and finds
% pairwise correlations between neurons over all trials. The second uses
% sorted_trains.mat and finds 7 correlation matrices between neurons
% depending on the TNR intensity.

num_neurons = 16;

% Load neuron trains
neuron_trains = load('neuron_trains.mat');
neuron_trains = neuron_trains.neuron_trains;
% Convert to matrix
neuron_trains = cell2mat(neuron_trains);
neuron_trains = transpose(neuron_trains);
% Compute pairwise correlations
corr = corrcoef(neuron_trains);
save('pairwise_corr.mat', 'corr');
% Set diagonal elements to 0
corr = corr - eye(num_neurons);
% Visualize correlation coefficients
image(corr, 'CDataMapping', 'scaled');
title('Pairwise Correlations');
colorbar;

% Load sorted trains
sorted_trains = load('sorted_trains.mat');
sorted_trains = sorted_trains.sorted_trains;

% Preallocate array for pairwise correlations across 7 TNR intensities
tnrs = [0; 60; 65; 70; 75; 80; 85];
num_tnrs = numel(tnrs);
tnr_corr = zeros(num_tnrs, num_neurons, num_neurons);

figure;
% For each TNR intensity...
for i = 1:num_tnrs
    % Grab trains for that intensity
    tnrtrains = sorted_trains(:,i);
    % Convert to matrix
    tnrtrains = cell2mat(tnrtrains);
    tnrtrains = transpose(tnrtrains);
    % Compute pairwise correlations
    corr = corrcoef(tnrtrains);
    tnr_corr(i,:,:) = corr;
    % Set diagonal elements to 0
    corr = corr - eye(num_neurons);
    % Visualize correlation coefficients
    subplot(3,3,i);
    image(corr, 'CDataMapping', 'scaled');
    titletext = ['Pairwise Correlations - TNR = ', num2str(tnrs(i))];
    title(titletext);
    % Print i 
    i
end

save('tnr_corr.mat', 'tnr_corr');
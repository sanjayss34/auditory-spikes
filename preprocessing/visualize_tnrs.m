% This file helps you visualize the TNRs across the 464 trials in order to
% better understand how to concatenate trials with similar TNRs.

% Load CohenNeurons.mat
CohenNeurons = load('CohenNeurons.mat');
CohenNeurons = CohenNeurons.CohenNeurons;

% Preallocate a 1x464 vector to store TNRs
TNRs = zeros(1,464);

% For each trial in Neuron 1 (identical across all neurons)...
for i = 1:num_trials
    if isnan(CohenNeurons(1).trials(i).TNR)
        TNRs(i) = 0;
    else
        TNRs(i) = CohenNeurons(1).trials(i).TNR;
    end
end

% Plot the TNRs to visualize the TNRs across trials
scatter(1:464, TNRs, 20, 'filled');
xlabel('Trial');
ylabel('TNR');
title('TNR vs. Trial');

% Plot histogram of TNRs
figure;
binranges = [0 57.5 62.5 67.5 72.5 77.5 82.5 87.5];
% Display histcounts
[histcounts, ind] = histc(TNRs, binranges);
% Make a table of TNR counts
counts1 = [0; 60; 65; 70; 75; 80; 85];
counts2 = transpose(histcounts(1:end-1))/464*100;
scatter(counts1, counts2, 20, 'filled');
xlabel('TNR');
ylabel('% of trials');
title('% of trials vs. TNR');

% ind indexes the trials by TNR
save('TNR_index.mat', 'ind');
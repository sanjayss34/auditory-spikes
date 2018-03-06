function [] = spike_train_generator(trial_length, bin_length)

% Create spike trains for all 16 x 464 trials
spike_train_processing(trial_length, bin_length);

% Sort trials by TNR and visualize the results
visualize_tnrs();

% Concatenate spike trains to get 1 train for every neuron at every TNR
concatenate_trains();

% Concatenate spike trains to get 1 train for every neuron
create_neuron_trains();

end
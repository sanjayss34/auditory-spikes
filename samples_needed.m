function [] = samples_needed()

    load 'neuron_trains.mat' neuron_trains;
    neuron_trains = cell2mat(neuron_trains);
    % Set sample sizes
    sample_sizes = [];
    % Set number of iterations to run estimate_ising
    iters = 1000;
    for i = 1:numel(sample_sizes)
        estimate_ising_core(iters, neuron_trains, sample_sizes(i));
    end

end
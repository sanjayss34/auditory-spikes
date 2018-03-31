function [] = plot_num_firing(h0, J, test_logical)
    load('neuron_trains.mat');
    neuron_trains = cell2mat(neuron_trains);
    neuron_trains = neuron_trains(:,test_logical);
    neuron_trains = (neuron_trains+1)/2;
    [N T] = size(neuron_trains);
    [sigm, states] = sample_ising_exact(h0, J);
    sigm = (sigm+1)/2;
    num_firing_each_sig = sum(sigm, 2);
    num_firing_bins = zeros(1, N+1);
    for i=1:numel(states)
        num_firing_bins(num_firing_each_sig(i)+1) = num_firing_bins(num_firing_each_sig(i)+1)+states(i);
    end
    h0_independent = log(mean(neuron_trains, 2)./(1-mean(neuron_trains, 2)))*0.5;
    h0_independent = transpose(h0_independent);
    [sigm, states] = sample_ising_exact(h0_independent, zeros(N, N));
    sigm = (sigm+1)/2;
    independent_num_firing_each_sig = sum(sigm, 2);
    independent_num_firing_bins = zeros(1, N+1);
    for i=1:numel(states)
        independent_num_firing_bins(independent_num_firing_each_sig(i)+1) = independent_num_firing_bins(independent_num_firing_each_sig(i)+1)+states(i);
    end
    experimental_num_firing = zeros(1, N+1);
    for i=1:size(neuron_trains, 2)
        experimental_num_firing(sum(neuron_trains(:, i))+1) = experimental_num_firing(sum(neuron_trains(:, i))+1)+1;
    end
    experimental_num_firing = experimental_num_firing/size(neuron_trains, 2);
    figure;
    plot(0:16, num_firing_bins);
    hold on
    plot(0:16, independent_num_firing_bins);
    plot(0:16, experimental_num_firing);
    hold off
    title('Probability vs. Number of neurons firing');
    xlabel('Number of neurons firing');
    ylabel('Probability');
    legend('Predicted (Ising)', 'Predicted (Independent)', 'Experimental');
end
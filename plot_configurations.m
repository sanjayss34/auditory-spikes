function [] = plot_configurations(h0, J, test_logical)
    load('neuron_trains');
    N = numel(h0);
    neuron_trains = cell2mat(neuron_trains);
    neuron_trains = neuron_trains(:,test_logical);
    neuron_trains2 = (neuron_trains+1)/2;
    [sigm, states] = sample_ising_exact(h0, J);
    h0_independent = log(mean(neuron_trains2, 2)./(1-mean(neuron_trains2, 2)))*0.5;
    h0_independent = transpose(h0_independent);
    [sigm, states_independent] = sample_ising_exact(h0_independent, zeros(N, N));
    experimental_probs = zeros(1, numel(states));
    for i=1:size(neuron_trains, 2)
        index = bi2de(transpose((neuron_trains(:,i)+1)/2))+1;
        experimental_probs(index) = experimental_probs(index)+1;
    end
    experimental_probs = experimental_probs/size(neuron_trains, 2);
    figure;
    scatter(experimental_probs, states);
    hold on;
    scatter(experimental_probs, states_independent);
    plot(linspace(10^(-5), 1, 1000), linspace(10^(-5), 1, 1000));
    set(gca, 'xscale', 'log');
    set(gca, 'yscale', 'log');
    title('Predicted vs. Experimental Probabilities for Configurations');
    xlabel('Experimental Probability');
    ylabel('Predicted Probability');
    legend('Ising Model', 'Independent', 'Identity function');
    hold off;
end
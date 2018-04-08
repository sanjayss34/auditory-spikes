function [h0 J train_logical] = estimate_ising_exact(iters)
    load 'neuron_trains.mat' neuron_trains;
    neuron_trains = cell2mat(neuron_trains);
    [N T] = size(neuron_trains);
    % Proportion used for training data
    p_train = 0.8;
    train_logical = false(T, 1);
    train_logical(1:round(p_train*T)) = true;
    train_logical = train_logical(randperm(T));
    test_logical = ~train_logical;
    neuron_trains_test = neuron_trains(:,test_logical);
    neuron_trains = neuron_trains(:,train_logical);
    
    mean_experiment = transpose(mean(neuron_trains,2));
    mean_experiment_product = neuron_trains*transpose(neuron_trains)/size(neuron_trains,2);

    % h0 = unifrnd(-1, 1, 1, N);
    neuron_trains2 = (neuron_trains+1)/2;
    % h0 = log(mean(neuron_trains2, 2)./(1-mean(neuron_trains2, 2)))*0.5;
    h0 = transpose(atanh(mean_experiment));
    h0 = transpose(h0);
    J = unifrnd(-0.1, 0.1, N, N);
    % J = zeros(N, N);
    maxdiff = 1;
    eta = 0.1;
    alpha = 0;
    prev_change_h0 = zeros(1, N);
    prev_change_J = zeros(N, N);
    itercount = 0;
    
    % Measure deviations from experiment
    sigma_diff = zeros(1,iters);
    corr_diff = zeros(1,iters);
    
    % Gradient Ascent
    sigm0 = 2*(randi(2,1,N)-1)-1;
    while itercount < iters
        
        % tic;
        % eta = eta/itercount;
        itercount = itercount+1;
        disp([itercount maxdiff]);
        maxdiff = 0;
        
        % Sample Ising Estimations
        sample_size = 1000;
        % [sigm, states] = sample_ising(sample_size, h0, J);
        % [sigm, states] = sample_ising_exact(h0, J);
        [sigm, states] = sw_sample_ising(h0, J, sample_size, sigm0);
        % if mod(itercount, 2) == 1
        
        % [sigm, states] = mh_sample_ising(1, sample_size, h0, J, 10, sigm0);
        size(sigm)
        size(states)
        sigm0 = sigm(sample_size,:);
        % end
        % [sigm, states] = gibbs_sample_ising(sample_size, h0, J, 100);
        weighted_states = sigm.*repmat(transpose(states), 1, size(sigm, 2));
        % toc
        
        % tic;
        % Update h0 
        mean_sigma = sum(weighted_states);
        diff = eta*(mean_experiment-mean_sigma) + alpha*prev_change_h0;
        prev_change_h0 = diff;
        h0 = h0 + diff;
        sigma_diff(itercount) = mean(abs(mean_experiment-mean_sigma));
        % toc
        
        % tic;
        % Update Jij
        mean_product = transpose(sigm)*weighted_states;
        diff = 0.5*eta*(mean_experiment_product-mean_product)+alpha*prev_change_J;
        diff(logical(eye(size(diff)))) = 0;
        maxdiff = max(max(max(maxdiff, abs(diff))));
        prev_change_J = diff;
        J = J + diff;
        corr_diff(itercount) = sum(sum(abs(mean_experiment_product-mean_product)))/(N^2);
        % toc
    end
    
    % Plot deviation from experiment over time
    figure;
    subplot(2,1,1);
    plot(1:itercount, sigma_diff(1:itercount));
    title('Deviation of Mean Firing Rate');
    subplot(2,1,2);
    plot(1:itercount, corr_diff(1:itercount));
    xlabel('# of Iterations');
    title('Deviation of Mean Correlation');

    mean_experiment = transpose(mean(neuron_trains_test,2));
    mean_experiment_product = neuron_trains_test*transpose(neuron_trains_test)/size(neuron_trains_test,2);    

    % Mean responses
    mrs = mean_sigma;
    mers = mean_experiment;
    % Mean products
    num_entries = N*(N-1)/2;
    meps = zeros(1,num_entries);
    mps = zeros(1,num_entries);
    k = 1;
    for i = 1:N
        for j = i+1:N
            meps(k) = mean_experiment_product(i,j) - mean_experiment(i)*mean_experiment(j);
            mps(k) = mean_product(i,j) - mean_sigma(i)*mean_sigma(j);
            k = k+1;
        end
    end
    
    % Plot predicted vs. empirical values
    figure;
    scatter(mers, mrs, 10, 'b', 'filled');
    hold on;
    xlabel('Mean Experimental Response');
    ylabel('Mean Predicted Response');
    lin = linspace(min(min(mers),min(mrs)),max(max(mers),max(mrs)),101);
    plot(lin, lin, 'r');
    figure;
    scatter(meps, mps, 10, 'b', 'filled');
    hold on;
    xlabel('Mean Experimental Correlation');
    ylabel('Mean Predicted Correlation');
    lin = linspace(min(min(meps),min(mps)),max(max(meps),max(mps)),101);
    plot(lin, lin, 'r');
end
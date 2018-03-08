function [h0 J] = estimate_ising(iters)
    load 'neuron_trains.mat' neuron_trains;
    neuron_trains = cell2mat(neuron_trains);
    [N T] = size(neuron_trains);
    h0 = unifrnd(-1, 1, 1, N);
    J = unifrnd(-1, 1, N, N);
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
    mean_experiment = transpose(mean(neuron_trains,2));
    mean_experiment_product = neuron_trains*transpose(neuron_trains)/size(neuron_trains,2);
    while itercount < iters
        
        tic;
        % eta = eta/itercount;
        itercount = itercount+1;
        disp([itercount maxdiff]);
        maxdiff = 0;
        
        % Sample Ising Estimations
        sample_size = 500;
        %samples = sample_ising(sample_size, h0, J);
        samples = mh_sample_ising(1, sample_size, h0, J, 500);
        toc
        
        tic;
        % Update h0 
        mean_sigma = mean(samples);
        diff = eta*(mean_experiment-mean_sigma) + alpha*prev_change_h0;
        prev_change_h0 = diff;
        h0 = h0 + diff;
        sigma_diff(itercount) = mean(abs(mean_experiment-mean_sigma));
        toc
        
        tic;
        % Update Jij
        mean_product = transpose(samples)*samples/size(samples,1);
        diff = eta*(mean_experiment_product-mean_product)+alpha*prev_change_J;
        diff(logical(eye(size(diff)))) = 0;
        maxdiff = max(max(max(maxdiff, abs(diff))));
        prev_change_J = diff;
        J = J + diff;
        corr_diff(itercount) = sum(sum(abs(mean_experiment_product-mean_product)))/(N^2);
        toc
    end
    
    % Plot deviation from experiment over time
    figure;
    subplot(2,1,1);
    plot(1:iters, sigma_diff);
    title('Deviation of Mean Firing Rate');
    subplot(2,1,2);
    plot(1:iters, corr_diff);
    xlabel('# of Iterations');
    title('Deviation of Mean Correlation');
    
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
            meps(k) = mean_experiment_product(i,j);
            mps(k) = mean_product(i,j);
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
    plot(lin, lin, 'r')
    figure;
    scatter(meps, mps, 10, 'b', 'filled');
    hold on;
    xlabel('Mean Experimental Correlation');
    ylabel('Mean Predicted Correlation');
    lin = linspace(min(min(meps),min(mps)),max(max(meps),max(mps)),101);
    plot(lin, lin, 'r')
end
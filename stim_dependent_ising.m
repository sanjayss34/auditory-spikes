function [h0 J] = stim_dependent_ising(iters)
    load 'sorted_trains.mat' sorted_trains;
    load 'neuron_trains.mat' neuron_trains;
    load('amplitude_trains.mat');
    % [N K] = size(sorted_trains);
    neuron_trains = cell2mat(neuron_trains);
    [N T] = size(neuron_trains);
    K = 10;
    binmarks = linspace(min(amplitude_trains), max(amplitude_trains), K+1);
    counts = zeros(1, K);
    mean_experiment = zeros(K, N);
    for t=1:T
        for j=1:K
            if amplitude_trains(t) >= binmarks(j) && amplitude_trains(t) < binmarks(j+1)
                counts(j) = counts(j)+1;
                mean_experiment(j,:) = mean_experiment(j,:)+transpose(neuron_trains(:,t));
            end
        end
    end
    for j=1:K
        mean_experiment(j,:) = mean_experiment(j,:)/counts(j);
    end
    mean_overall_experiment = mean(transpose(neuron_trains));
    h0 = repmat(atanh(mean_overall_experiment), K, 1);
    J = unifrnd(-1, 1, N, N);
    maxdiff = 1;
    eta = 0.1;
    alpha = 0;
    prev_change_h0 = zeros(K, N);
    prev_change_J = zeros(N, N);
    itercount = 0;

    % Measure deviations from experiment
    sigma_diff = zeros(1,iters);
    corr_diff = zeros(1,iters);

    % Gradient Ascent
    mean_experiment_product = neuron_trains*transpose(neuron_trains)/size(neuron_trains,2);
    while itercount < iters

        % eta = eta/itercount;
        itercount = itercount+1;
        disp([itercount maxdiff]);
        maxdiff = 0;
        
        % Sample Ising Estimations
        sample_size = 1000;
        mean_product = zeros(N, N);
        mean_sigma_all = zeros(K, N);
        for k=1:K
            % tic;
            [sigm, states] = sample_ising_exact(h0(k,:), J);
            weighted_states = sigm.*transpose(states);
            % toc
        
            % tic;
            % Update h0 
            mean_sigma = sum(weighted_states);
            mean_sigma_all(k,:) = mean_sigma;
            diff = eta*(mean_experiment(k,:)-mean_sigma) + alpha*prev_change_h0(k,:);
            maxdiff = max(max(abs(diff)), maxdiff);
            prev_change_h0(k,:) = diff;
            h0(k,:) = h0(k,:) + diff;
            sigma_diff(itercount) = mean(abs(mean_experiment(k,:)-mean_sigma));
            % toc
            mean_product = mean_product+counts(k)*(transpose(sigm)*weighted_states)/sum(counts);
        end
        
        tic;
        % Update Jij
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
    plot(1:itercount, sigma_diff);
    title('Deviation of Mean Firing Rate');
    subplot(2,1,2);
    plot(1:itercount, corr_diff);
    xlabel('# of Iterations');
    title('Deviation of Mean Correlation');
    
    % Mean responses
    mrs = mean_sigma_all(K,:);
    mers = mean_experiment(K,:);
    
    overall_mean_experiment = mean(neuron_trains, 2);
    overall_mean_model = zeros(N, 1);
    for k=1:K
        overall_mean_model = overall_mean_model+mean_sigma_all(k,:)*counts(k)/sum(counts);
    end
    % Mean products
    num_entries = N*(N-1)/2;
    meps = zeros(1,num_entries);
    mps = zeros(1,num_entries);
    k = 1;
    for i = 1:N
        for j = i+1:N
            meps(k) = mean_experiment_product(i,j) - overall_mean_experiment(i)*overall_mean_experiment(j);
            mps(k) = mean_product(i,j) - overall_mean_model(i)*overall_mean_model(j);
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
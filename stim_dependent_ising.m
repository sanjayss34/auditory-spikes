function [h0 J train_logical] = stim_dependent_ising(iters)
    load 'sorted_trains.mat' sorted_trains;
    load 'neuron_trains.mat' neuron_trains;
    load('feature_intensity.mat');
    % amplitude_trains = cell2mat(feature_intensity);
    % [N K] = size(sorted_trains);
    % neuron_trains = cell2mat(neuron_trains);
    amplitude_trains = feature_intensity;
    [N T] = size(neuron_trains);
    % Proportion used for training data
    p_train = 1-50/T;
    train_logical = false(T, 1);
    train_logical(1:round(p_train*T)) = true;
    train_logical = train_logical(randperm(T));
    test_logical = ~train_logical;
    neuron_trains_test = neuron_trains(:,test_logical);
    neuron_trains = neuron_trains(:,train_logical);
    amplitude_trains_test = amplitude_trains(:,test_logical);
    amplitude_trains = amplitude_trains(:,train_logical);
    K = 10;
    binmarks = linspace(min(min(amplitude_trains)), max(max(amplitude_trains)), K+1);
    counts = zeros(1, K);
    mean_experiment = zeros(K, N);
    for n=1:N
        for k=1:K
            mean_experiment(k,n) = mean(neuron_trains(n,(amplitude_trains(n,:) >= binmarks(k)) & (amplitude_trains(n,:) < binmarks(k+1))));
        end
    end
    time_bin_map = zeros(T, N);
    for n=1:N
        for k=1:K
            time_bin_map(amplitude_trains(n,:) >= binmarks(k) & amplitude_trains(n,:) < binmarks(k+1),n) = k;
        end
    end
%     for t=1:T
%         for j=1:K
%             mean_experiment(j,
%             for n=1:N
%                 if amplitude_trains(n, t) >= binmarks(j) && amplitude_trains(n, t) < binmarks(j+1)
%                     counts(j) = counts(j)+1;
%                     mean_experiment(j,n) = mean_experiment(j,n)+transpose(neuron_trains(:,t));
%                 end
%             end
%         end
%     end
%     for j=1:K
%         mean_experiment(j,:) = mean_experiment(j,:)/counts(j);
%     end
    mean_overall_experiment = mean(transpose(neuron_trains));
    h0 = transpose(repmat(atanh(mean_overall_experiment), K, 1));
    J = unifrnd(-1, 1, N, N);
    maxdiff = 1;
    eta = 0.1;
    alpha = 0;
    prev_change_h0 = zeros(N, K);
    prev_change_J = zeros(N, N);
    itercount = 0;

    % Measure deviations from experiment
    sigma_diff = zeros(1,iters);
    corr_diff = zeros(1,iters);

    % Gradient Ascent
    mean_experiment_product = neuron_trains*transpose(neuron_trains)/size(neuron_trains,2);
    time_points = randi([1 sum(train_logical)], 1, iters);
    train_indices = find(train_logical);
    while itercount < iters

        % eta = eta/itercount;
        itercount = itercount+1;
        % disp([itercount maxdiff]);
        maxdiff = 0;
        
        t = train_indices(time_points(itercount));
        h0_t = zeros(1, N);
        for n=1:N
            h0_t(n) = h0(n,time_bin_map(t,n));
        end
        h0_t
        J
        [sigm, states] = sample_ising_exact(h0_t, J);
        weighted_states = sigm.*transpose(states);
        mean_sigma = sum(weighted_states);
        diff = eta*(neuron_trains(:,t)-mean_sigma);
        maxdiff = max(max(abs(diff)), maxdiff);
        for n=1:N
            h0(n,time_bin_map(t,n)) = h0(n,time_bin_map(t,n))+diff(n);
        end
        pred_product = transpose(sigm)*weighted_states;
        actual_product = neuron_trains(:,t)*transpose(neuron_trains(:,t));
        diff = eta*(actual_product-pred_product)+alpha*prev_change_J;
        diff(logical(eye(size(diff)))) = 0;
        maxdiff = max(maxdiff, max(max(abs(diff))));
        prev_change_J = diff;
        J = J + diff;
        
%         % Sample Ising Estimations
%         sample_size = 1000;
%         mean_product = zeros(N, N);
%         mean_sigma_all = zeros(K, N);
%         for k=1:K
%             % tic;
%             [sigm, states] = sample_ising_exact(h0(k,:), J);
%             weighted_states = sigm.*transpose(states);
%             % toc
%         
%             % tic;
%             % Update h0 
%             mean_sigma = sum(weighted_states);
%             mean_sigma_all(k,:) = mean_sigma;
%             diff = eta*(mean_experiment(k,:)-mean_sigma) + alpha*prev_change_h0(k,:);
%             maxdiff = max(max(abs(diff)), maxdiff);
%             prev_change_h0(k,:) = diff;
%             h0(k,:) = h0(k,:) + diff;
%             sigma_diff(itercount) = mean(abs(mean_experiment(k,:)-mean_sigma));
%             % toc
%             mean_product = mean_product+counts(k)*(transpose(sigm)*weighted_states)/sum(counts);
%         end
%         
%         tic;
%         % Update Jij
%         diff = eta*(mean_experiment_product-mean_product)+alpha*prev_change_J;
%         diff(logical(eye(size(diff)))) = 0;
%         maxdiff = max(max(max(maxdiff, abs(diff))));
%         prev_change_J = diff;
%         J = J + diff;
%         corr_diff(itercount) = sum(sum(abs(mean_experiment_product-mean_product)))/(N^2);
%         toc
    end
    
    % Plot deviation from experiment over time
%     figure;
%     subplot(2,1,1);
%     plot(1:itercount, sigma_diff);
%     title('Deviation of Mean Firing Rate');
%     subplot(2,1,2);
%     plot(1:itercount, corr_diff);
%     xlabel('# of Iterations');
%     title('Deviation of Mean Correlation');
    
    % Mean responses
    mean_experiment = transpose(mean(neuron_trains_test,2));
    mean_experiment_product = neuron_trains_test*transpose(neuron_trains_test)/size(neuron_trains_test,2);
    mean_model = zeros(1, N);
    mean_product_model = zeros(N, N);
    for n=1:N
        for k=1:K
            time_bin_map(amplitude_trains_test(n,:) >= binmarks(k) & amplitude_trains_test(n,:) < binmarks(k+1)) = k;
        end
    end
    for t=1:numel(test_logical)
        if test_logical(t) == 1
            h0_t = zeros(1, N);
            for n=1:N
                h0_t(n) = h0(n,time_bin_map(t,n));
            end
            [sigm, states] = sample_ising_exact(h0_t, J);
            weighted_states = sigm.*transpose(states);
            mean_sigma = sum(weighted_states);
            pred_product = transpose(sigm)*weighted_states;
            mean_model = mean_model+mean_sigma;
            mean_product_model = mean_product_model+pred_product;
        end
    end
    mean_model = mean_model/sum(test_logical);
    mean_product_model = mean_product_model/sum(test_logical);
    mrs = mean_model;
    mers = mean_experiment;

%     mrs = mean_sigma_all(K,:);
%     mers = mean_experiment(K,:);
    
%     overall_mean_experiment = mean(neuron_trains, 2);
%     overall_mean_model = zeros(N, 1);
%     for k=1:K
%         overall_mean_model = overall_mean_model+mean_sigma_all(k,:)*counts(k)/sum(counts);
%     end
    % Mean products
    num_entries = N*(N-1)/2;
    meps = zeros(1,num_entries);
    mps = zeros(1,num_entries);
    k = 1;
    for i = 1:N
        for j = i+1:N
            meps(k) = mean_experiment_product(i,j) - mean_experiment(i)*mean_experiment(j);
            mps(k) = mean_product_model(i,j) - mean_model(i)*mean_model(j);
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
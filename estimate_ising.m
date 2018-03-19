function [h0 J] = estimate_ising()
    load('neuron_trains.mat');
    neuron_trains = cell2mat(neuron_trains);
    [N T] = size(neuron_trains);
    experimental_means = zeros(N, 1);
    corr = zeros(N, N);
    for i=1:N
        experimental_means(i) = mean(neuron_trains(i,:));
    end
    load('sorted_trains.mat');
    [N, K] = size(sorted_trains);
    experimental_means_stim = zeros(N, K);
    for j=1:K
        for i=1:N
            experimental_means_stim(i, j) = mean(cell2mat(sorted_trains(i, j)));
        end
    end
    for i=1:N
        for j=1:N
            if i == j
                corr(i, j) = 1;
            else
                corr(i, j) = mean(neuron_trains(i,:).*neuron_trains(j,:))-experimental_means(i)*experimental_means(j);
            end
        end
    end
    h0 = unifrnd(0, 1, 1, N);
    J = unifrnd(0, 1, N, N);
    maxdiff = 1;
    eta = 0.1;
    alpha = 0;
    prev_change_h0 = zeros(1, N);
    prev_change_J = zeros(N, N);
    itercount = 1;
    avg_mean_deviations = [];
    avg_corr_deviations = [];
    min_max_diff = maxdiff;
    best_h0 = h0;
    best_J = J;
    while itercount < 51
        if mod(itercount, 50) == 0
            eta = eta/2;
        end
        disp([itercount maxdiff]);
        itercount = itercount+1;
        maxdiff = 0;
        prev_h0 = h0;
        prev_J = J;
        sample_size = 500;
        % samples = mh_sample_ising(1, ones(1, N), sample_size, h0, J, 500);
        samples = sample_ising(sample_size, h0, J);
        avg_mean_deviation = 0;
        for i=1:N
            mean_sigma_i = mean(samples(:,i));
            deviation = experimental_means(i)-mean_sigma_i;
            diff = eta*(deviation)+alpha*prev_change_h0(i);
            avg_mean_deviation = avg_mean_deviation+abs(deviation);
            maxdiff = max(maxdiff, abs(deviation));
            prev_change_h0(i) = diff;
            h0(i) = h0(i)+diff;
        end
        avg_mean_deviation = avg_mean_deviation/N;
        avg_corr_deviation = 0;
        for i=1:N
            for j=i+1:N
                mean_experiment_product = corr(i, j)+experimental_means(i)*experimental_means(j);
                mean_product = mean(samples(:,i).*samples(:,j));
                deviation = mean_experiment_product-mean_product;
                diff = eta*(deviation)+alpha*prev_change_J(i,j);
                avg_corr_deviation = avg_corr_deviation+abs(deviation-(experimental_means(i)*experimental_means(j)-mean(samples(:,i))*mean(samples(:,j))));
                maxdiff = max(maxdiff, abs(deviation));
                prev_change_J(i, j) = diff;
                prev_change_J(j, i) = diff;
                J(i, j) = J(i,j)+diff;
                J(j, i) = J(j,i)+diff;
            end
        end
        avg_corr_deviation = avg_corr_deviation/(N*(N+1)/2);
        avg_mean_deviations = [avg_mean_deviations avg_mean_deviation];
        avg_corr_deviations = [avg_corr_deviations avg_corr_deviation];
        if maxdiff < min_max_diff
            min_max_diff = maxdiff;
            best_h0 = prev_h0;
            best_J = prev_J;
        end
    end
    h0 = best_h0;
    J = best_J;
    subplot(2, 1, 1);
    plot(5:numel(avg_mean_deviations), avg_mean_deviations(5:numel(avg_mean_deviations)));
    xlabel('Iteration');
    ylabel('Avg. deviation');
    title('Avg. deviation of mean response vs. Iteration number');
    subplot(2, 1, 2);
    plot(1:numel(avg_corr_deviations), avg_corr_deviations);
    xlabel('Iteration');
    ylabel('Avg. deviation');
    title('Avg. deviation of correlation vs. Iteration number');
    figure;
    samples = sample_ising(sample_size, h0, J);
    pred_means = [];
    for i=1:N
        pred_mean = mean(samples(:,i));
        pred_means = [pred_means pred_mean];
    end
    scatter(experimental_means, pred_means);
    hold on;
    plot(linspace(-1, 0, 101), linspace(-1, 0, 101));
    xlabel('Experimental mean');
    ylabel('Predicted mean');
    title('Predicted vs. Experimental means');
    hold off;
    figure;
    pred_corrs = [];
    exp_corrs = [];
    for i=1:N
        for j=i+1:N
            pred_corr = mean(samples(:,i).*samples(:,j))-mean(samples(:,i))*mean(samples(:,j));
            pred_corrs = [pred_corrs pred_corr];
            exp_corrs = [exp_corrs corr(i, j)];
        end
    end
    scatter(exp_corrs, pred_corrs);
    hold on;
    plot(linspace(-0.05, 0.15, 101), linspace(-0.05, 0.15, 101));
    xlabel('Experimental correlation');
    ylabel('Predicted correlation');
    title('Predicted vs. Experimental correlations');
    hold off;
end
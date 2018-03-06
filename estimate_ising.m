function [h0 J] = estimate_ising()
    load('neuron_trains.mat');
    neuron_trains = cell2mat(neuron_trains);
    [N T] = size(neuron_trains);
    h0 = unifrnd(0, 1, 1, N);
    J = unifrnd(0, 1, N, N);
    maxdiff = 1;
    eta = 1;
    alpha = 0;
    prev_change_h0 = zeros(1, N);
    prev_change_J = zeros(N, N);
    itercount = 1;
    while maxdiff > 0.01
        eta = eta/sqrt(itercount);
        disp([itercount maxdiff]);
        itercount = itercount+1;
        maxdiff = 0;
        sample_size = 500;
        samples = sample_ising(sample_size, h0, J);
        for i=1:N
            mean_sigma_i = mean(samples(:,i));
            mean_experiment_i = mean(neuron_trains(i,:));
            diff = eta*(mean_experiment_i-mean_sigma_i)+alpha*prev_change_h0(i);
            maxdiff = max(maxdiff, abs(diff));
            prev_change_h0(i) = diff;
            h0(i) = h0(i)+diff;
        end
        for i=1:N
            for j=i+1:N
                mean_experiment_product = mean(neuron_trains(i,:).*neuron_trains(j,:));
                mean_product = mean(samples(:,i).*samples(:,j));
                diff = eta*(mean_experiment_product-mean_product)+alpha*prev_change_J(i,j);
                maxdiff = max(maxdiff, abs(diff));
                prev_change_J(i, j) = diff;
                prev_change_J(j, i) = diff;
                J(i, j) = J(i,j)+diff;
                J(j, i) = J(j,i)+diff;
            end
        end
    end
end
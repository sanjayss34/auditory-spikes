function [std_h, std_J, std_corr, chi] = standard_deviations(h0, J)
load('neuron_trains.mat');
neuron_trains = cell2mat(neuron_trains);
[N T] = size(neuron_trains);
chi = zeros(N*(N+1)/2, N*(N+1)/2);
[sigm, states] = sample_ising_exact(h0, J);
weighted_states = sigm.*repmat(transpose(states), 1, size(sigm, 2));
means = sum(weighted_states);
mean_products = transpose(sigm)*weighted_states;
% mean_triple_products = zeros(N, N, N);
% for i=1:size(sigm, 1)
%     disp(sprintf("A %d", i));
%     for j1=1:N
%         for j2=1:N
%             for j3=1:N
%                 mean_triple_products(j1, j2, j3) = mean_triple_products(j1, j2, j3) + states(i)*sigm(i,j1)*sigm(i,j2)*sigm(i,j3);
%             end
%         end
%     end
% end
% mean_quartuple_products = zeros(N, N, N, N);
% for i=1:size(sigm, 1)
%     disp(sprintf("B %d", i));
%     for j1=1:N
%         for j2=1:N
%             for j3=1:N
%                 for j4=1:N
%                     mean_quartuple_products(j1, j2, j3, j4) =mean_quartuple_products(j1, j2, j3, j4)+states(i)*sigm(i,j1)*sigm(i,j2)*sigm(i,j3)*sigm(i,j4);
%                 end
%             end
%         end
%     end
% end
% save('mean_triple_products.mat', "mean_triple_products");
% save('mean_quartuple_products.mat', "mean_quartuple_products");
load('mean_triple_products.mat');
load('mean_quartuple_products.mat');
for i=1:N
    for j=1:N
        chi(i, j) = mean_products(i, j)-means(i)*means(j);
    end
end
for i=1:N
    count = N+1;
    for j=1:N
        for k=1:(j-1)
            chi(i, count) = mean_triple_products(i, j, k)-means(i)*mean_products(j, k);
            chi(count,i) = chi(i, count);
            count=count+1;
        end
    end
end
count1 = N+1;
for i=1:N
    for j=1:(i-1)
        count2 = N+1;
        for k=1:N
            for l=1:(k-1)
                chi(count1, count2) = mean_quartuple_products(i, j, k, l)-mean_products(i, j)*mean_products(k, l);
                chi(count2, count1) = chi(count1, count2);
                count2=count2+1;
            end
        end
        count1=count1+1;
    end
end
inv_chi = inv(chi);
std_h = zeros(1, N);
for i=1:N
    std_h(i) = sqrt(inv_chi(i, i)/size(neuron_trains, 2));
end
std_J = zeros(N, N);
count = N+1;
for i=1:N
    for j=1:(i-1)
        std_J(i, j) = sqrt(inv_chi(count, count)/size(neuron_trains, 2));
        std_J(j, i) = std_J(i, j);
        count=count+1;
    end
end
std_corr = zeros(N, N);
for i=1:N
    for j=1:(i-1)
        std_corr(i, j) = sqrt((mean_products(i, j)*(1-mean_products(i, j)))/size(neuron_trains, 2));
        std_corr(j, i) = std_corr(i, j);
    end
end
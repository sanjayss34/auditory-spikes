function [sigm, states2]  = sample_ising(m, h0, J)
    N = numel(h0);
    states = zeros(1, 2^N);
    i = 1:numel(states);
    sigm = de2bi(i-1,N);
    sigm = sigm*2-1;
    for i=1:numel(states)
        exponential_arg = sum(h0.*sigm(i,:))+0.5*sum(sum(J.*(transpose(sigm(i,:))*sigm(i,:))));
        states(i) = exp(exponential_arg);
    end
    states = states/sum(states);
    distribution = makedist('Multinomial', 'probabilities', states);
    %sigm = zeros(m, N);
    %{
    for j=1:m
        sample = random(distribution);
        sigm(j,:) = de2bi(sample-1, N);
    end
    %}
    sample = random(distribution, m, 1);
    sigm = de2bi(sample-1, N);
    sigm = 2*sigm-1;
    %states2 = ones(1,m)/m;
    states2 = states(sample);
    states2 = states2/sum(states2);
end
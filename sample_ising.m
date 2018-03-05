function [samples]  = sample_ising(m, h0, J)
    N = numel(h0);
    states = zeros(1, 2^N);
    for i=1:numel(states)
        sigm = de2bi(i-1, N);
        sigm = sigm*2-1;
        exponential_arg = sum(h0.*sigm)+0.5*sum(sum(J.*(transpose(sigm)*sigm)));
        states(i) = exp(exponential_arg);
    end
    states = states/sum(states);
    distribution = makedist('Multinomial', 'probabilities', states);
    samples = zeros(m, N);
    for j=1:m
        sample = random(distribution);
        samples(j,:) = de2bi(sample-1, N);
    end
    samples = 2*samples-1;
end
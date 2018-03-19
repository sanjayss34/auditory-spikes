function [sigm, states]  = sample_ising_exact(h0, J)
    %h0 is a row vector
    %J is a square matrix
    %sigm is a 2^N-x-N matrix
    N = numel(h0);
    states = zeros(1, 2^N);
    numbers = transpose(linspace(1,numel(states),numel(states)));
    sigm = de2bi(numbers-1)*2-1;
    exponential_arg = zeros(1, 2^N);
    % Zero out the diagonal of J
    J(logical(eye(size(J)))) = 0;
    for i=1:numel(states)
        exponential_arg(i) = sum(h0.*sigm(i,:))+0.5*sum(sum(J.*(transpose(sigm(i,:))*sigm(i,:))));
    end
    states = exp(exponential_arg);
    states = states/sum(states);
end
function [sigm, states] = single_chain_mh2(beta, m, h0, J, burn_iters, interval_iters, sigm0)  
    N = numel(h0);    
    % If sigma is not specified, assume all neurons are +1
    if ~exist('sigm0', 'var')
        sigm0 = unifrnd(-1,1,1,N);
    end
    iters = burn_iters+m*interval_iters;
    selected_neurons = randi([1 N], 1, iters);
    
    sigm = zeros(iters,N);
    for t = 1:iters
        % Zero out the diagonal of J
        J(logical(eye(size(J)))) = 0;
        % Current energy
        exponential_arg_u = -1*(sum(h0.*sigm0)-sum(sum(J.*(transpose(sigm0)*sigm0))));
        Hu = exponential_arg_u;
        % Flipped energy
        temp_sigm = sigm0;
        temp_sigm(selected_neurons(t)) = -temp_sigm(selected_neurons(t));
        exponential_arg_v = -1*(sum(h0.*temp_sigm)-sum(sum(J.*(transpose(temp_sigm)*temp_sigm))));
        Hv = exponential_arg_v;
        % Compare energies
        if Hv > Hu
            A = exp(-beta*(Hv-Hu));
            if rand < A
                sigm0 = temp_sigm;
            end
        else
            sigm0 = temp_sigm;
        end
        sigm(t,:) = sigm0;
    end
    
    ind = (1:m).*interval_iters+burn_iters;
    sigm = sigm(ind,:);
    states = ones(1,m)/m;
end
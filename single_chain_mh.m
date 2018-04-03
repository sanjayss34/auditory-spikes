function [sigm, states] = single_chain_mh(beta, m, h0, J, burn_iters, interval_iters, sigm0)  
    N = numel(h0);
    sigm = zeros(m, N);    
    
    % If sigma is not specified, assume all neurons are +1
    if ~exist('sigm0', 'var')
        sigm0 = unifrnd(-1,1,1,N);
    end
    iters = burn_iters+m*interval_iters;
    selected_neurons = randi([1 N], 1, iters);
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
        % Save sample at every interval after burn-in period
        if mod(t-burn_iters, interval_iters) == 0 && t > burn_iters
            ind = (t-burn_iters)/interval_iters;
            sigm(ind,:) = sigm0;
        end
    end
    states = ones(1,m)/m;
end
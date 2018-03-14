function [samples] = single_chain_mh(beta, m, h0, J, burn_iters, interval_iters, sigm0)  
    N = numel(h0);
    samples = zeros(m, N);    
    
    % If sigma is not specified, assume all neurons are +1
    if ~exist('sigm0', 'var')
        sigm0 = unifrnd(-1,1,1,N);
    end
    sigm = sigm0;
    selected_neurons = randi([1 N], 1, iters);
    iters = burn_iters+m*interval_iters;
    for t = 1:iters
        % Current energy
        exponential_arg_u = -1*(sum(h0.*sigm)+0.5*sum(sum(J.*(transpose(sigm)*sigm))));
        Hu = exponential_arg_u;
        % Flipped energy
        temp_sigm = sigm;
        temp_sigm(selected_neurons(t)) = -temp_sigm(selected_neurons(t));
        exponential_arg_v = -1*(sum(h0.*temp_sigm)+0.5*sum(sum(J.*(transpose(temp_sigm)*temp_sigm))));
        Hv = exponential_arg_v;
        % Compare energies
        if (Hv-Hu) > 0
            A = exp(-beta*(Hv-Hu));
            if rand < A
                sigm = temp_sigm;
            end
        else
            sigm = temp_sigm;
        end
        % Save sample at every interval after burn-in period
        if mod(t-burn_iters, interval_iters) == 0
            ind = (t-burn_iters)/interval_iters;
            samples(ind,:) = sigm;
        end
    end
end
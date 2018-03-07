function [samples]  = mh_sample_ising(beta, sigm0, m, h0, J, iters)  
    N = numel(h0);
    samples = zeros(m, N);    
    for i = 1:m
        % If sigma is not specified, assume all neurons are +1
        if ~exist('sigm0', 'var')
            sigm0 = ones(1,N);
        end
        sigm = sigm0;
        selected_neurons = randi([1 N], 1, iters);
        for t = 1:iters
            % Current energy
            exponential_arg = sum(h0.*sigm)+0.5*sum(sum(J.*(transpose(sigm)*sigm)));
            Hu = exp(exponential_arg);
            % Flipped energy
            temp_sigm = sigm;
            temp_sigm(selected_neurons(t)) = -temp_sigm(selected_neurons(t));
            exponential_arg = sum(h0.*temp_sigm)+0.5*sum(sum(J.*(transpose(temp_sigm)*temp_sigm)));
            Hv = exp(exponential_arg);
            % Compare energies
            if (Hv-Hu) > 0
                A = exp(-beta*(Hv-Hu));
                if rand < A
                    sigm = temp_sigm;
                end
            else
                sigm = temp_sigm;
            end
        end
        samples(i,:) = sigm;
    end
end
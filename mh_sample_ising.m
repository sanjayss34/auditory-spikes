function [sigm, states] = mh_sample_ising(beta, m, h0, J, iters)  
    N = numel(h0);
    sigm = zeros(m, N);    
    for i = 1:m
        % If sigma is not specified, assume all neurons are +1
        if ~exist('sigm0', 'var')
            sigm0 = 2*(randi(2,1,N)-1)-1;
        end
        sigm_i = sigm0;
        selected_neurons = randi([1 N], 1, iters);
        for t = 1:iters
            % Current energy
            exponential_arg_u = -1*(sum(h0.*sigm_i)+0.5*sum(sum(J.*(transpose(sigm_i)*sigm_i))));
            Hu = exponential_arg_u;
            % Flipped energy
            temp_sigm = sigm_i;
            temp_sigm(selected_neurons(t)) = -temp_sigm(selected_neurons(t));
            exponential_arg_v = -1*(sum(h0.*temp_sigm)+0.5*sum(sum(J.*(transpose(temp_sigm)*temp_sigm))));
            Hv = exponential_arg_v;
            % Compare energies
            if (Hv-Hu) > 0
                A = exp(-beta*(Hv-Hu));
                if rand < A
                    sigm_i = temp_sigm;
                end
            else
                sigm_i = temp_sigm;
            end
        end
        sigm(i,:) = sigm_i;
    end
    
    states = ones(1,m)/m;
end
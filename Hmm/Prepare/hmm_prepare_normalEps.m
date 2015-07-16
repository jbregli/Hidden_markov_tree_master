function [epsilon] = hmm_prepare_normalEps(epsilon)
% hmm_prepare_normalEps: Normalize 'epsilon' so that each row sums to 1.

    %% Initilisation:
    d_eps = length(size(epsilon));

    %% Epsilon:
    if d_eps == 2     
        n_state = size(epsilon,1);
        
        % Normalize:
        epsilon = epsilon ./ repmat(sum(epsilon,2),1,n_state);   

    else
        n_state = size(epsilon,3);
        
        % Normalize:
        for row=1:n_state
            tmp_sum = sum(epsilon(:,:,row,:),4);
            for col=1:n_state
                epsilon(:,:,row,col) =epsilon(:,:,row,col) ./ tmp_sum;
            end
        end
    end
end


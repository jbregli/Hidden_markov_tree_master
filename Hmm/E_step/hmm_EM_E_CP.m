function [S] = hmm_EM_E_CP(S, theta)
% hmm_EM_E_CP: COMPUTE THE DESIRED CONDITIONAL PROBABILITIES
%   
%   See "Wavelet-based statistical signal processing using hidden markov
%   models" for more details on the algorithm
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'hmm_prepare'.
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%
%   --------
%   OUTPUTS:
%   --------
%   - S: cell(struct)
%       Structure updated
%   -pmf_state: (?)
%       Probability mass function for the hidden state variables.
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Lot of optimization possible in the computation of the betas
%   (beta_exclude node is computed several times, useless for loops (?),
%   useless variables...

    %% Initialization:
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);
    end
        
    %% Conditional probabilities:
    % Loop over all layers:
    for layer=1:n_layer
        % Loop over all scales:
        for scale=1:n_elmt(1,layer)
            temp_sum = hmm_cp_sum(S, layer, scale);                         % ----------> Check OK
            S = hmm_cp_givenVis(S, temp_sum, layer, scale);                 % ----------> Check OK
            S = hmm_cp_givenVisAndParents(S, theta, temp_sum, layer, scale);
        end
    end
end


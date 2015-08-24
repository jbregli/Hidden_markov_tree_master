function [set_S] = hmm_EM_E(set_S, theta)
% hmm_EM_E: COMPUTE THE EXPECTATION STEP OF THE "EXPECTATION/MAXIMISATION"
%           ALGORITHM
%   
%   See "Wavelet-based statistical signal processing using hidden markov
%   models" for more details on the algorithm
%
%   --------
%   INPUTS:
%   --------
%   - set_S: cell{cell(struct)}
%       Set of structures obtained with the function 'hmm_EM_E'.      
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%
%   --------
%   OUTPUTS:
%   --------
%   - set_S: cell(struct)
%       Set of structures updated
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Lot of optimization possible in the computation of the betas
%   (beta_exclude node is computed several times, useless for loops (?),
%   useless variables...

    %% Initialization:
    
    %% Upward-Downward algorithm
    % Loop over the set of scattering transform:
    for im=1:length(set_S)       
        % Up step: Transmits information from last layer to the first and
        % computes the betas.
        set_S{im} = hmm_EM_E_up(set_S{im}, theta);                          % ----------> Check OK
             
        % Down step: Transmits information from first layer to the last and
        % computes the alphas
        set_S{im} = hmm_EM_E_down(set_S{im}, theta);                        % ----------> Check OK  
        
        % Compute the conditional probabilities:
        set_S{im} = hmm_EM_E_CP(set_S{im}, theta);                          % ----------> Check OK   
        
    end
    
end


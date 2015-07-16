function [S] = hmm_EM_E_up(S, theta)
% hmm_EM_E_upward: UP PASS OF THE UPWARD/DOWNWARD ALGORITHM
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
%   --------s
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
    
    % For all the states variable s_i of leaf nodes,
    % calculate Beta_i(m) for each state m:
    for layer = n_layer:-1:1
        for scale=1:n_elmt(1,layer)
            if isempty(S{layer}.hmm{scale}.children)
                % Beta_i(m) - S{n_layer}.hmm{scale}.beta.excludeNode :
                S = hmm_beta_givenNode(S, theta, layer, scale);
            end
        end
    end
    
    %% Loop over the layers (bottom up loop)
    for layer= n_layer:-1:2
%         Beta_i(m): - orignal algo out of the loop - some nodes have no
%         children and thus it may need to be inside the loop
%         for scale=1:n_elmt(1,layer)
%             S = hmm_beta_givenNode(S, theta, layer, scale);                 % ----------> Check OK
%         end
        
        % beta_{i,rho(i)}(m) - S{layer}.hmm{scale}.beta.givenParents :
        for scale=1:n_elmt(1,layer)
            S = hmm_beta_givenParents(S, theta, layer, scale);              % ----------> Check Durand and Crouse have a different formula here.
                                                                            %             Durand version is more convincing.
        end
        
        % beta_{rho(i)}(m) - S{layer}.hmm{scale}.beta.Parents :
        for scale=1:n_elmt(1,layer)
            S = hmm_beta_Parents(S, theta, layer, scale);                   % ----------> Check OK
                                                                            %             ERROR HERE -> the result of this should be use to update 
                                                                            %             S{layer-1}.hmm{'parents(scale)'}.beta.givenNode            
        end
        
        % beta_{rho(i)\i}(m) - S{layer}.hmm{scale}.beta.excludeNode
        for scale=1:n_elmt(1,layer)
            S = hmm_beta_excludeNode(S, theta, layer, scale);               % ----------> Check OK
        end
        
    end
end


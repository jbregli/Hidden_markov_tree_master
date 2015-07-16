function [S] = hmm_EM_E_down(S, theta)
% hmm_EM_E_upward: DOWN STEP OF THE UPWARD/DOWNWARD ALGORITHM
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
    for l=1:n_layer
        n_elmt(1,l) = length(S{l}.signal);
    end
    n_state = size(theta{1}.proba{1}, 3);
        
    %% Loop over the layers:
    % First layer:  
    S{1}.hmm{1}.alpha = theta{1}.proba{1};
        
    for layer=2:n_layer
        % alpha_i(m) = sum_{n=1}^{M} (epsilon_{i,rho(i)}^{mn} alpha_{rho(i)}(n) beta_{rho(i)\i}(n) ):
        for scale=1:n_elmt(1,layer)
            % Reset to 0 to be sure:
            S{layer}.hmm{scale}.alpha = zeros(size(S{layer}.hmm{scale}.alpha));
            
            % All possible state for this node:
            for m=1:n_state
                % Sum over all possible states for the father:
                for n=1:n_state
                    S{layer}.hmm{scale}.alpha(:,:,m) = ...
                        S{layer}.hmm{scale}.alpha(:,:,m) + ...
                        theta{layer}.epsilon{scale}(:,:, n, m) ...
                            .* S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.alpha(:,:,n) ...
                            .* S{layer}.hmm{scale}.beta.excludeNode(:,:,n);
                end
            end
            
%             % +++ Correction of elements falsly rounded-up to 0 - avoid overflow:
%             S{layer}.hmm{scale}.alpha(S{layer}.hmm{scale}.alpha==0) = 0.0001;    
%             
            %% +++ Sanity check:
            if max(max(max(isnan( S{layer}.hmm{scale}.alpha))))
                disp(['Hmm_Down: NAN alpha at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
            if max(max(max(S{layer}.hmm{scale}.alpha ==0))) ==1
                disp(['alpha: 0 alpha at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
            % +++            
            
        end
    end
    
    %% +++ Sanity check:
    if max(max(max(isnan( S{1}.hmm{1}.alpha))))
        disp(['Hmm_Down: NAN alpha at layer ' num2str(1) ...
            ' and scale ' num2str(1)])
    end
    if max(max(max(S{1}.hmm{1}.alpha ==0))) ==1
        disp(['Hmm_Down: 0 alpha at layer ' num2str(1) ...
            ' and scale ' num2str(1)])
    end
    % +++
    
end


function [S] = hmm_cp_givenVisAndParents(S, theta, sum_norm, layer, scale)
% hmm_beta: COMPUTE THE CONDITIONAL PROBABILITY OF A STATE GIVEN THE 
%           SCATTERING COEFFICIENTS.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'hmm_prepare'. 
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%   - layer: int
%       Layer to be considered.
%   - scale: int
%       Scale to be considered.

%   --------
%   OUTPUTS:
%   --------
%   - temp: float
%       Normalization constant
%   --------
%   IMPROVEMENTS:
%   --------

    %% Initialization:
    % Sizes:
    n_state = size(theta{1}.proba{1}, 3);   

    %% Compute the conditional probability given the scattering coefficient:
    % Not define for first layer --> =1
    if layer > 1
        % Loop over child's states:
        for m=1:n_state
            % Loop over father's states:
            for n=1:n_state
%                 disp(['+++ hmm-cp_givenVisAndParents: sub2ind(n,m) = ' num2str( sub2ind([n_state,n_state], n, m))])
%                 disp(['+++ hmm-cp_givenVisAndParents: m = ' num2str(m)])
%                 disp(['+++ hmm-cp_givenVisAndParents: n = ' num2str(n)])
                
                S{layer}.hmm{scale}.condProb.givenVisAndParents(:,:, n, m) = ...
                          S{layer}.hmm{scale}.beta.givenNode(:,:,m) ...
                            .* theta{layer}.epsilon{scale}(:,:, n, m) ...
                            .* S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.alpha(:,:,n)...
                            .* S{layer}.hmm{scale}.beta.excludeNode(:,:,n) ...
                                ./ sum_norm;
            end
        end

        % +++ Correction of elements falsly rounded-up to 0 - avoid overflow:
        S{layer}.hmm{scale}.condProb.givenVisAndParents(S{layer}.hmm{scale}.condProb.givenVisAndParents==0) = 0.0001;      
        
        %% +++ Sanity check:
        if max(max(max(isnan(S{layer}.hmm{scale}.condProb.givenVisAndParents))))
            disp(['CP - givenVisAndParents: NAN at layer ' num2str(layer) ...
                ' and scale ' num2str(scale)])
        end
        if max(max(max(S{layer}.hmm{scale}.condProb.givenVisAndParents == 0))) == 1
            disp(['CP - givenVisAndParents: 0 at layer ' num2str(layer) ...
                ' and scale ' num2str(scale)])
        end
        % +++
    end

end


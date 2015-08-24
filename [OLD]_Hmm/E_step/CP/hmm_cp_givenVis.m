function [S] = hmm_cp_givenVis(S, sum_norm, layer, scale)
% hmm_beta: COMPUTE THE CONDITIONAL PROBABILITY OF A STATE GIVEN THE 
%           SCATTERING COEFFICIENTS.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'hmm_prepare'.  
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
    n_state = size(S{1}.hmm{1}.alpha, 3);
    
    % Create a temp_sum variable of size = [size(sum_norm) n_state] where
    % temp_sum(:,:,i) = sum_norm
    temp_sum = zeros([size(sum_norm) n_state]);
    for m=1:n_state
        temp_sum(:,:,m) = sum_norm;
    end
    
    %% Compute the conditional probability given the scattering coefficient:
    S{layer}.hmm{scale}.condProb.givenVis = ...
            S{layer}.hmm{scale}.alpha .* S{layer}.hmm{scale}.beta.givenNode ...
                ./ temp_sum;

    % +++ Correction of elements falsly rounded-up to 0 - avoid overflow:
    S{layer}.hmm{scale}.condProb.givenVis(S{layer}.hmm{scale}.condProb.givenVis==0) = 0.0001;      
          
    %% +++ Sanity check:      
    if max(max(max(isnan(S{layer}.hmm{scale}.condProb.givenVis))))
        disp(['CP - givenVis: NAN at layer ' num2str(layer) ...
              ' and scale ' num2str(scale)])
    end
    if max(max(max(S{layer}.hmm{scale}.condProb.givenVis == 0))) ==1
        disp(['CP - givenVis: 0 at layer ' num2str(layer) ...
              ' and scale ' num2str(scale)])
    end
    % +++   

end


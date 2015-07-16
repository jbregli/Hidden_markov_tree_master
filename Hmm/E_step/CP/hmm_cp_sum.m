function [sum_norm] = hmm_cp_sum(S, layer, scale)
% hmm_beta: COMPUTE THE NORMALIZATION CONSTANT FOR THE CONDITIONAL 
%           PROBABILITIES  
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
    
    %% Compute the Normalisation sum:
    temp_prod = S{layer}.hmm{scale}.alpha .* S{layer}.hmm{scale}.beta.givenNode;
    sum_norm = sum(temp_prod, 3);

%     % +++ Correction of elements falsly rounded-up to 0 - avoid overflow:
%     sum_norm(sum_norm==0) = 0.0001;  
                                            
%     % +++ Sanity check:
%     if max(max(isnan(sum_norm)))
%         disp(['cp_sum: NAN sum_norm at layer ' num2str(layer) ...
%               ' and scale ' num2str(scale)])
%     end
%     if max(max( sum_norm == 0)) ==1
%         disp(['cp_sum: 0 sum_norm at layer ' num2str(layer) ...
%               ' and scale ' num2str(scale)])
%     end
%     % +++   
end


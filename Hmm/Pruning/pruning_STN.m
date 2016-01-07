function [signal_to_noise] = pruning_STN(theta)

% pruning_HMT_STN: PRUNE THE TREE ACCORDING TO SIGNAL TO NOISE RATION.
%
%   THis function compute the signal to noise ratio  of each leaf of the
%   tree and remove 'prun_param.rmv_prc'% of the leaf with the highest 
%   signal to noise ratio (\frac{\abs(\mu)}{\sigma}) .
%
%   --------
%   INPUTS:
%   --------
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib used
%       to store the model parameters.
%   - prun_params: (optional) struct
%       Each fields is a meta parameter for the pruning.
%           - .rmv_prc: (optional) float(0,1) (default: 0.5)
%               Percentage of highest signal to noise ratio leafs to be
%               removed
%           - .n_iteration: (optional) int (default: 1)
%               Number of iteration over the pruning. If all the children 
%               of a node are removed during the first pruning iteration
%               then the second pass is done on the   ---- TBC
%       
%   --------
%   OUTPUTS:
%   --------
%   - theta: cell(struct)
%       Structure with a new organisation as some node are now ignored.
%
%   --------
%   TODO:
%   --------

    %% Preparation:          
    % Sizes:
    n_layer = length(theta);
    n_scale = zeros(1,n_layer);
    
    signal_to_noise.tree = cell(1,n_layer);
    signal_to_noise.list = [];
        
    for layer=1:n_layer
        % Number of scale: per layer:
        n_scale(1,layer) = length(theta{layer}.proba);
        
        signal_to_noise.tree{layer} = cell(1,n_scale(1,layer));
    end
    
    %% Signal to noise ratio:
    for layer=1:n_layer
        for scale=1:n_scale(1,layer);
            tmp_stn = min(theta
            signal_to_noise.tree
            
        end
    end
end


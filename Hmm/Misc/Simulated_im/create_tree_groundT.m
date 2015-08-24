function [S] = create_tree_groundT(S, theta_gt)
% GT_theta:  CREATE THE STRUCTURE TO STORE THE TREE TRUE PARAMETERS.
%   Given a scattering transform (cell of structure), this functions 
%   associates to each node the required fields for the HMM modeling.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Set of structures obtained with the function 'scat' of the 
%       'scatnet' lib.     
%   - theta_gt: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%       theta{layer}.proba{index}
%       theta{layer}.epsilon{index}
%       theta{layer}.mu{index}
%       theta{layer}.sigma{index}
%       theta{layer}.distr{index}
%
%   --------
%   OUTPUTS:
%   --------
%   - theta_gt: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%       theta{layer}.proba{index}
%       theta{layer}.epsilon{index}
%       theta{layer}.mu{index}
%       theta{layer}.sigma{index}
%       theta{layer}.distr{index}
%
%   --------
%   IMPROVEMENTS:
%   --------
%

    %% Initialization:   
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);
    end
           
    %% Sample a tree given the parameters in theta:
    for layer=1:n_layer
        for scale=1:n_elmt(1,layer)
            % Sample the state:
            S = sample_state(S, theta_gt, layer, scale);
            % Sample 'image value':
            S = sample_image(S, theta_gt, layer, scale);   
        end
    end
            
            
    
    

end


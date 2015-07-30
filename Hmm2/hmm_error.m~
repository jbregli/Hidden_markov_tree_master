function error_LW = hmm_error(theta_est, theta_gt) %, cv_ach_strct)
% hmm_error: COMPUTE THE SQUARE ERROR BETWEEN THE ESTIMATE AND THE GROUND
%            TRUTH
%   Given the estimated theta and the ground truth, this function compute
%   the error. If a 'cv_achieved_strct' is provided, this function computes
%   2 errors. One over all the pixels and the other over all the pixels
%   that have already converged.
%
%   --------
%   INPUTS:
%   --------
%   - theta_est: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib
%       containing the estimated modelisation parameters.
%   - theta_gt: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib
%       containing the ground truth modelisation parameters.

%   --------
%   OUTPUTS:
%   --------
%   - error: Struct
%       Fields: .mu , .sigma , .epsilon, . proba
%       Error over all the layers
%   - error_LW: cell(Struct)
%       Fields: .mu , .sigma , .epsilon, . proba
%       Note that some fields will be unused at some layer (only 'proba' is
%       used at layer one, while it's not used for any other layer)
%
%   --------
%   ISSUES:
%   --------
%
%   --------
%   TODO:
%   --------

    %% Preparation:
    % optional variables:

    % Sizes:
    n_layer = length(theta_est);
    s_image = size(theta_est{1}.proba{1}(:,:,1));
    n_pixel = prod(s_image);
    n_state = size(theta_est{1}.proba{1},3);

    n_scale = zeros(1,n_layer);
    
    % Structure to store errors layer wise:
    error_LW = cell(1, n_layer);
    
    for layer=1:n_layer
        n_scale(1,layer) = length(theta_est{layer}.proba);

        % Structure:
        error_LW{layer}.proba = cell(1,n_scale(1,layer));
        error_LW{layer}.mu = cell(1,n_scale(1,layer));
        error_LW{layer}.sigma = cell(1,n_scale(1,layer));           
        error_LW{layer}.epsilon= cell(1,n_scale(1,layer));
    end
    
    %% Error:
    % Layer 1 - proba:
    tmp_proba_N = (theta_est{1}.proba{1} - theta_gt{1}.proba{1}).^2;
    tmp_proba_N = sum(sum(sum(tmp_proba_N,3),2),1);
    
    tmp_proba_D = (theta_est{1}.proba{1} + theta_gt{1}.proba{1}).^2;
    tmp_proba_D = sum(sum(sum(tmp_proba_D,3),2),1);
    
    error_LW{1}.proba{1} = tmp_proba_N / tmp_proba_D;
    
    % Loop over the layers:
    for layer=2:n_layer 
        %  Loop over the scales at 'layer':
        for scale=1:n_scale(1,layer)
            % Mu:
            tmp_mu = (theta_est{layer}.mu{scale} - ...
                theta_gt{layer}.mu{scale}).^2;
            tmp_mu = sum(sum(sum(tmp_mu,3),2),1);
            error_LW{layer}.mu{scale} = tmp_mu / n_pixel;
            
            % Sigma:
            tmp_sigma = (theta_est{layer}.sigma{scale} - ...
                theta_gt{layer}.sigma{scale}).^2;
            tmp_sigma = sum(sum(sum(tmp_sigma,3),2),1);
            error_LW{layer}.sigma{scale} = tmp_sigma / n_pixel;
            
            % Epsilon:
            tmp_eps = (theta_est{layer}.epsilon{scale} - ...
                theta_gt{layer}.epsilon{scale}).^2;
            tmp_eps = sum(sum(sum(sum(tmp_eps,4),3),2),1);
            error_LW{layer}.epsilon{scale} = tmp_eps / (2*n_pixel);
        end
    end
end


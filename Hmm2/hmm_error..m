function [cv_ach_strct, cv_ach_bool] = ...
    hmm_error(theta_est, theta_gt)
% hmm_error: COMPUTE THE ERROR BETWEEN THE ESTIMQTE AND THE GROUND TRUTH


    %% Preparation:
    % optional variables:
    % add cv achieved structure to compute error overall and error over
    % converged pixels
    
    % Sizes:
    n_layer = length(theta_est);
    n_scale = zeros(1,n_layer);
    s_image = size(theta_est{1}.proba{1}(:,:,1));
    n_state = size(theta_est{1}.proba{1},3);

    % Structure to store the 'proba'
    error = cell(1, n_layer);
    
    for layer=1:n_layer
        n_scale(1,layer) = length(theta{layer}.proba);

        % Structure:
        error{layer}.proba = cell(1,n_scale(1,layer));
        error{layer}.mu = cell(1,n_scale(1,layer));
        error{layer}.sigma = cell(1,n_scale(1,layer));
        error{layer}.epsilon = cell(1,n_scale(1,layer));
    end
    
    %% Error:
    % Loop over the layers:
    for layer=1:n_layer 
        % Loop over the scales at 'layer':
        for scale=1:n_scale(1,layer)
            % COMPUTE THE ERROR (SQUARE ERROR)
            % ADD A SECOND ERROR FOR CONVERGED PIXELS
        end
    end
end


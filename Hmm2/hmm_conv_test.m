function [cv_ach_struct, cv_ach_bool] = ...
    hmm_conv_test(theta, theta_old, step, mixing, sensibility)
% hmm_conv_test: CONVERGENCE TEST FOR EM ALGORITHM
%   Given a theta at step n and theta at step n-1, this function asseses if
%   convergence has occcured.

    %% Preparation:
    % optional variables:
    if ~exist('mixing','var')
        mixing = 10;
    end
    if ~exist('sensibility','var')
        sensibility = 1e-4;
    end

    % Sizes:
    n_layer = length(theta);
    n_elmt = zeros(1,n_layer);
    s_image = size(theta{1}.proba{1}(:,:,1));
    n_state = size(theta{1}.proba{1},3);

    % Structure to store convergence status:
    cv_ach_struct = cell(1, n_layer);
    
    for layer=1:n_layer
        n_elmt(1,layer) = length(theta{layer}.proba);

        % Structure:
        cv_ach_struct{layer}.proba = cell(1,n_elmt(1,layer));
        cv_ach_struct{layer}.epsilon = cell(1,n_elmt(1,layer));
        cv_ach_struct{layer}.mu = cell(1,n_elmt(1,layer));
        cv_ach_struct{layer}.sigma = cell(1,n_elmt(1,layer));
        
        % Initialization of the matrices:
        for i=1:n_elmt(1,layer)
            cv_ach_struct{layer}.proba{i} = zeros([s_image n_state]);
            cv_ach_struct{layer}.epsilon{i} = zeros([s_image n_state n_state]); 
            cv_ach_struct{layer}.mu{i} = zeros([s_image n_state]);
            cv_ach_struct{layer}.sigma{i} = zeros([s_image n_state]);
        end
    end
    
    cv_ach_bool = true;
    
    fields = fieldnames(cv_ach_struct{1});

    %% Convergence test:
    % If the number of step is too small then cv didn't occur (burning
    % time)
    if step < mixing
        cv_ach_bool = false;
    else
        % Loop over the layers:
        for layer=1:n_layer
            %  Loop over the scales at 'layer':
            for scale=1:n_elmt(1,layer)
                for i=1:numel(fields)
                    % Delta:
                    cv_ach_struct{layer}.(fields{i}){scale} = ...
                        theta{layer}.(fields{i}){scale} ...
                        - theta_old{layer}.(fields{i}){scale};
                    
                    % Convergence structure boolean:
                    cv_ach_struct{layer}.(fields{i}){scale}(...
                        cv_ach_struct{layer}.(fields{i}){scale} < sensibility) = true;
                    cv_ach_struct{layer}.(fields{i}){scale}(...
                        cv_ach_struct{layer}.(fields{i}){scale} >= sensibility) = false;
                    
                    % Convergence boolean:
                    if any(any(any(cv_ach_struct{layer}.(fields{i}){scale}) == 0))
                        cv_ach_bool = false;
                    end                        
                end
            end
        end
    end
end


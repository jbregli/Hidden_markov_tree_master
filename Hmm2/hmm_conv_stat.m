function cv_stat = hmm_conv_stat(cv_ach_strct)
% hmm_conv_stat: CONVERGENCE TEST FOR EM ALGORITHM
%   Given a theta at step n and theta at step n-1, this function asseses if
%   convergence has occcured.

cv_lim = 5; % TB ADDED AS INPUT

    %% Preparation:
    % optional variables:
    
    % Sizes:
    n_layer = length(cv_ach_strct);
    s_image = size(cv_ach_strct{1}.proba{1}(:,:,1));
    n_pixel = prod(s_image);
    n_state = size(cv_ach_strct{1}.proba{1},3);
    n_scale = zeros(1,n_layer);
    

    % Structure to store convergence error:
    cv_stat = cell(1, n_layer);

    for layer=1:n_layer   
        if layer ==  1
            n_scale(1,layer) = length(cv_ach_strct{layer}.proba);
            
            % Structure:
            cv_stat{layer}.proba = cell(1,n_scale(1,layer));
            cv_stat{layer}.global = 0;

        else
            n_scale(1,layer) = length(cv_ach_strct{layer}.mu);
            
            % Structure:
            cv_stat{layer}.epsilon = cell(1,n_scale(1,layer));
            cv_stat{layer}.mu = cell(1,n_scale(1,layer));
            cv_stat{layer}.sigma = cell(1,n_scale(1,layer));
            cv_stat{layer}.global = 0;

        end
    end
       
    %% Convergence error:
    % Loop over the layers:
    for layer=1:n_layer 
        fields = fieldnames(cv_stat{layer});    
        
        %  Loop over the scales at 'layer':
        for scale=1:n_scale(1,layer)
            for i=1:numel(fields)                     
                if ~strcmp(fields{i},'global') 
                    % Count:
                    tmp_count = ...
                        cv_ach_strct{layer}.(fields{i}){scale} > cv_lim;

                    % Error count:
                    if strcmp(fields{i},'epsilon')
                        cv_stat{layer}.(fields{i}){scale} = ...
                            sum(sum(sum(sum(tmp_count,4),3),2),1) ...
                            / n_pixel;
                    else
                        cv_stat{layer}.(fields{i}){scale} = ...
                            sum(sum(sum(tmp_count,3),2),1) / n_pixel;
                    end
                           
                    % Global error:
                    cv_stat{layer}.global = cv_stat{layer}.global ...
                        + cv_stat{layer}.(fields{i}){scale};                                 
                end
            end
        end
        cv_stat{layer}.global = cv_stat{layer}.global ...
            / ((numel(fields) -1) * n_scale(1,layer));
    end
end


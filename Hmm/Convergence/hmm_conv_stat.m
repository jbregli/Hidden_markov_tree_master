function cv_stat = hmm_conv_stat(cv_status)
% hmm_conv_stat: CONVERGENCE TEST FOR EM ALGORITHM
%   Given a theta at step n and theta at step n-1, this function asseses if
%   convergence has occcured.


    %% Preparation:
    % optional variables:
    
    % Sizes:
    n_layer = length(cv_status.params);
    s_image = size(cv_status.params{1}.proba{1}(:,:,1));
    n_pixel = prod(s_image);
    n_state = size(cv_status.params{1}.proba{1},3);
    n_scale = zeros(1,n_layer);
    

    % Structure to store convergence error:
    cv_stat.overallCv = 0.;
    cv_stat.layerwise = cell(1, n_layer);

    for layer=1:n_layer   
        if layer ==  1
            n_scale(1,layer) = length(cv_status.params{layer}.proba);
            
            % Structure:
            cv_stat.layerwise{layer}.proba = cell(1,n_scale(1,layer));
            cv_stat.layerwise{layer}.global = 0;

        else
            n_scale(1,layer) = length(cv_status.params{layer}.mu);
            
            % Structure:
            cv_stat.layerwise{layer}.epsilon = cell(1,n_scale(1,layer));
            cv_stat.layerwise{layer}.mu = cell(1,n_scale(1,layer));
            cv_stat.layerwise{layer}.sigma = cell(1,n_scale(1,layer));
            cv_stat.layerwise{layer}.global = 0;

        end
    end
       
    %% Convergence error:
    % Loop over the layers:
    for layer=1:n_layer 
        fields = fieldnames(cv_stat.layerwise{layer});    
        
        %  Loop over the scales at 'layer':
        for scale=1:n_scale(1,layer)
            for i=1:numel(fields)                     
                if ~strcmp(fields{i},'global') 
                    % Count:
                    tmp_count = ...
                        cv_status.params{layer}.(fields{i}){scale} > cv_status.limit;

                    % Error count:
                    if strcmp(fields{i},'epsilon')
                        cv_stat.layerwise{layer}.(fields{i}){scale} = ...
                            sum(sum(sum(sum(tmp_count,4),3),2),1) ...
                            / n_pixel;
                    else
                        cv_stat.layerwise{layer}.(fields{i}){scale} = ...
                            sum(sum(sum(tmp_count,3),2),1) / n_pixel;
                    end
                           
                    % Global error:
                    cv_stat.layerwise{layer}.global = cv_stat.layerwise{layer}.global ...
                        + cv_stat.layerwise{layer}.(fields{i}){scale};                                 
                end
            end
        end
        cv_stat.layerwise{layer}.global = cv_stat.layerwise{layer}.global ...
            / ((numel(fields) -1) * n_scale(1,layer));
        
        % Overall conv:
        cv_stat.overallCv = cv_stat.overallCv ...
            + cv_stat.layerwise{layer}.global /n_layer;        
    end
end


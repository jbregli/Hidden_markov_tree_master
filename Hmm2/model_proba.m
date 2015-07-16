function [Pvalue] = model_proba(S, theta, layer, scale, threshold)
% model_proba: Summary of this function goes here
%   Detailed explanation goes here

    %% Preparation:
     % Sizes:
    n_state = size(theta{1}.proba{1}, 3);
    s_im = size(S{1}.signal{1});
    
    Pvalue = zeros([s_im n_state]);
    
    % Threshold:
    if nargin < 5
        b_threshold = false;
    else
        b_threshold = true;
    end
    
    %% Compute the PDF:
    if strcmp(theta{layer}.distr,'MixtGauss')        
        for m=1:n_state           
            Pvalue(:,:,m) = normpdf(...
                  S{layer}.signal{scale},...
                  squeeze(theta{layer}.mu{scale}(:,:,m)),...
                  squeeze(theta{layer}.sigma{scale}(:,:,m))); 
        end      
        
        % +++ Try to avoid overflow by tresholding:
        if b_threshold
            Pvalue = Pvalue.*(Pvalue>threshold) + threshold*(Pvalue<=threshold);
        end
    end           
end


function [S] = sample_image(S, theta_gt, layer, scale)
% sample_state: SAMPLE THE IMAGE VALUE OF THE SIMULATED ST.
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Set of structures obtained with the function 'scat' of the 
%       'scatnet' lib. 
%   - theta_gt: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%   - layer: int
%       Layer to be considered.
%   - scale: int
%       Scale to be considered.
%
%   --------
%   OUTPUTS:
%   --------
%   - S: cell(struct)
%       Structure updated
%
%   --------
%   IMPROVEMENTS:
%   --------
%

    %% Initialization:     
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1,n_layer);
    for l=1:n_layer
        n_elmt(1,l) = length(S{l}.signal);
    end
    n_state = size(theta_gt{1}.proba{1}, 3);  
    s_im = size(S{1}.signal{1});

    %% Sample image:
    mu = zeros(s_im);
    sigma = zeros(s_im);
    
    for m=1:n_state
        s_index = (S{layer}.hmm{scale}.state == m);

        tmp_mu = theta_gt{layer}.mu{scale}(:,:,m);
        tmp_sigma = theta_gt{layer}.sigma{scale}(:,:,m);
               
        mu(s_index) = tmp_mu(s_index);
        sigma(s_index) = tmp_sigma(s_index);
    end
        
    % Sample "pixel" values from a gaussian with the given params:
    S{layer}.signal{scale} = normrnd(mu, sigma);
end



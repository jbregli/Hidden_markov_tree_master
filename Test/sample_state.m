function [S] = sample_state(S, theta_gt, layer, scale)
% sample_state: SAMPLE THE STATE OF THE SIMULATED ST.
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

    %% Sample states:
    if layer == 1
        % Sampling from the root node state distribution:
        uni = rand(s_im);
        
        cumprob = zeros([s_im, (n_state+1)]);
        cumprob(:,:,2:end) = cumsum(theta_gt{layer}.proba{scale},3);
        
        
        for n=1:n_state
            ind= (uni>cumprob(:,:,n)) & (uni<=cumprob(:,:,n+1));
            S{layer}.hmm{scale}.state(ind)= n;
        end
        
    else
        % State is conditioned by the father node's state:
        f_layer = layer-1;
        f_scale = S{layer}.hmm{scale}.parent;
        f_state = S{f_layer}.hmm{f_scale}.state;
        
        % Loop over the fathers' states:
        for n=1:n_state
            idx = (f_state == n);
            
            trs_prob = theta_gt{layer}.epsilon{scale}(n,:);
            
            uni=rand(length(idx));
            cumprob = zeros([s_im, (n_state+1)]);
            cumprob(:,:,2:end) = cumsum(theta_gt{2}.epsilon{1}(:,:,n:n_state:end));
            
            % disp(['+++ sample_state: n = ' num2str(n) , ' epsilon(n) = ' num2str(theta_gt{2}.epsilon{1}(1,1,n:n_state:end))])
            
            % Loop over the sons' states
            for m = 1:n_state
                ind= (uni>cumprob(:,:,m)) & (uni<=cumprob(:,:,m+1));
                S{layer}.hmm{scale}.state(idx(ind)) = m;
            end
        end
    end
end



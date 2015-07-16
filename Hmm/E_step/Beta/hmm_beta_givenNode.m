function [S] = hmm_beta_givenNode(S, theta, layer, scale)
% hmm_beta: COMPUTE THE DENSITY FUNCTION OF EACH STATE WITH THE GIVEN 
%           PARAMETERS  
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'hmm_prepare'.
%   - theta: cell(struct)
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
%   - Add 'assert' for the distribution

    %% Initialization:
    % Sizes:
    n_state = size(theta{1}.proba{1}, 3);
    
    %% Compute the PDF:
    if strcmp(theta{layer}.distr,'MixtGauss')
        for m=1:n_state         
            S{layer}.hmm{scale}.beta.givenNode(:,:,m) = normpdf(...
                                  S{layer}.signal{scale},...
                                  squeeze(theta{layer}.mu{scale}(:,:,m)),...
                                  squeeze(theta{layer}.sigma{scale}(:,:,m))); 
        end       
    end                 

    % +++ Prevent overfolwing:
    S{layer}.hmm{scale}.beta.givenNode(abs(S{layer}.hmm{scale}.beta.givenNode) < 1e-4) = 1e-4;
    
    %% +++ Sanity check:
    if max(max(max(isnan( S{layer}.hmm{scale}.beta.givenNode))))
        disp(['Hmm_beta -givenNode : NAN at layer ' num2str(layer) ...
            ' and scale ' num2str(scale)])
    end
    if max(max(max(S{layer}.hmm{scale}.beta.givenNode ==0))) ==1
        disp(['Hmm_beta - givenNode : 0 at layer ' num2str(layer) ...
            ' and scale ' num2str(scale)])
        S{layer}.signal{scale}
        disp(['    mu= ' num2str(min(min(min(theta{layer}.mu{scale}))))])
        disp(['    sigma= ' num2str(min(min(min(theta{layer}.sigma{scale}))))])

    end
    % +++         
    
    
end


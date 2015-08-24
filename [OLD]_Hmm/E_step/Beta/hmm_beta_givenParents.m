function [S] = hmm_beta_givenParents(S, theta, layer, scale)
% hmm_beta: COMPUTE THE DENSITY FUNCTION OF EACH STATE GIVEN ITS PARENTS
%  
%   We have:
%   beta_{i,rho(i)}(m) = sum_{n=1}^{M}(epsilon_{i,rho(i)}^{mn} beta_i(n))
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

%   --------
%   OUTPUTS:
%   --------
%   - S: cell(struct)
%       Structure updated
%
%   --------
%   IMPROVEMENTS:
%   --------

    %% Initialization:
    % Sizes:
    n_state = size(theta{1}.proba{1}, 3);
    
    %% Compute the PDF:
    % beta_{i,rho(i)}(m) = sum_{n=1}^{M}(epsilon_{i,rho(i)}^{mn} beta_i(n))
    % Reset the matrix to 0:
    S{layer}.hmm{scale}.beta.givenParents = zeros(size(S{layer}.hmm{scale}.beta.givenParents));  
    
    for m=1:n_state                    
        % Sum over all possible states for the father:
        for n=1:n_state           
            S{layer}.hmm{scale}.beta.givenParents(:,:,m) = ...
                    S{layer}.hmm{scale}.beta.givenParents(:,:,m) ...
                    + ...
                    theta{layer}.epsilon{scale}(:,:,n, m) ...
                               .* S{layer}.hmm{scale}.beta.givenNode(:,:,n);                % WARNING: implementation of the version from "Durand".
                                                                                            % In epsilon this is the state of the children that change.
        end       
    end

%     % +++ Correction of elements falsly rounded-up to 0 - avoid overflow:
%     S{layer}.hmm{scale}.beta.givenParents(S{layer}.hmm{scale}.beta.givenParents==0) = 0.0001;    

    % +++ Sanity check:
    if max(max(max(isnan( S{layer}.hmm{scale}.beta.givenParents))))
        disp(['hmm_beta - givenParents: NAN at layer ' num2str(layer) ...
              ' and scale ' num2str(scale)])
    end
    if max(max(max(S{layer}.hmm{scale}.beta.givenParents == 0 ))) == 1
        disp(['hmm_beta - givenParents: 0 at layer ' num2str(layer) ...
              ' and scale ' num2str(scale)])
    end
    % +++
end


function [S] = hmm_beta_excludeNode(S, theta, layer, scale)
% hmm_beta: COMPUTE THE DENSITY FUNCTION OF THE SUBTREE STATE GIVEN ITS
%           PARENT
%
%   We have:
%   beta_{rho(i)\i}(m) = beta_{rho(i)}(m) / beta_{i,rho(i)}(m)
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
    
    
    %% Compute the PDF:
    % +++ ISSUE: Will generate NaN values if S{layer}.hmm{scale}.beta.givenParents
    % is equal to 0
    % Quick patch: add a small delta to givenParents to avoid the 0.
    % Better solution to be found   
    
    S{layer}.hmm{scale}.beta.excludeNode = ...
              S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.beta.givenNode ...
                    ./ S{layer}.hmm{scale}.beta.givenParents;                             

    %% +++ Sanity check:
    if max(max(max(isnan( S{layer}.hmm{scale}.beta.excludeNode))))
        disp(['Hmm_beta - excludeNode : NAN at layer ' num2str(layer) ...
            ' and scale ' num2str(scale)])
    end
    if max(max(max(S{layer}.hmm{scale}.beta.excludeNode ==0))) ==1
        disp(['Hmm_beta - excludeNode : 0 at layer ' num2str(layer) ...
            ' and scale ' num2str(scale)])
    end
    % +++                                
                
end


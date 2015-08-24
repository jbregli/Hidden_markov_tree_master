function [S] = hmm_beta_Parents(S, theta, layer, scale)
% hmm_beta: COMPUTE THE DENSITY FUNCTION OF EACH PARENT
%
%   We have:
%   beta_{rho(i)}(m) = g(w_{rho(i)}; mu_{rho(i),m}sigma_{rho(i),m}) 
%                             x prod_{i in c(rho(i))} beta_{i,rho(i)}(m)
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

    %% Initialization:
    
    %% Compute the PDF:
    % beta_{rho(i)}(m) = g(w_{rho(i)}; mu_{rho(i),m}sigma_{rho(i),m}) 
    %                            x prod_{i in c(rho(i))} beta_{i,rho(i)}(m)
    
    % Product of all beta_given_parents of the layer:
    tmp_prod = ones(size(S{layer}.hmm{scale}.beta.givenParents));
    
    % Product over all the brothers:
    for i=1:length(S{layer}.hmm{scale}.brothers)
        tmp_prod = tmp_prod .* ...
           S{layer}.hmm{S{layer}.hmm{scale}.brothers(i)}.beta.givenParents;
    end
    
    % Likelihood of the father node:
    S = hmm_beta_givenNode(S, theta, (layer-1), S{layer}.hmm{scale}.parent);
   
    % Update the value of 'beta.givenNode' of the parent:
    S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.beta.givenNode = ...
              S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.beta.givenNode ...
                  .* tmp_prod;
 
    %% +++ Sanity check:
    if max(max(max(isnan(S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.beta.givenNode))))
        disp(['Hmm_beta - parents : NAN at layer_c ' num2str(layer) ...
            ' and scale_c ' num2str(scale) ...
            ' layer_f ' num2str(layer-1) ...
            ' and scale_f ' num2str(S{layer}.hmm{scale}.parent)])
    end
    if max(max(max(S{layer-1}.hmm{S{layer}.hmm{scale}.parent}.beta.givenNode ==0))) ==1
        disp(['Hmm_beta - parents : 0 at layer_c ' num2str(layer) ...
            ' and scale_c ' num2str(scale) ...
            ' layer_f ' num2str(layer-1) ...
            ' and scale_f ' num2str(S{layer}.hmm{scale}.parent)])
    end
    % +++          
              
end


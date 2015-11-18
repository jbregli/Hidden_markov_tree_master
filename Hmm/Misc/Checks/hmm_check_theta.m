function [check_bool] = hmm_check_theta(theta)
%hmm_Scheck_0nan: PERFORM A SANITY CHECK ON THE GIVEN VARIABLE 
%   Is there any 0, NaN, or infinite values.
%
%   --------
%   INPUTS:
%   --------
%   - var_TBC: Multidimensional Arrays (d=3 or 4)
%       Matrix of size [s_im, n_state] or [s_im, n_state, n_state] holding 
%       the value to check.
%   - var_fnct: string
%       Function name from where the checking is done. Help for debugging.
%       Display only if 'verbose' is true.
%   - var_name: string
%       'Var_TBC' name. Help for debugging. Display only if 'verbose' is
%       true.
%   - layer: int
%       Layer number of 'var_TBC'. Help for debugging. Display only if 
%       'verbose' is true.
%   - scale: int
%       Scale number of 'var_TBC'. Help for debugging. Display only if 
%       'verbose' is true.
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%   - nan: (optional) bool (default= true)
%       Check for nan values?
%   - zero: (optional) bool (default= true)
%       Check for 0 values?
%   - infty: (optional) bool (default= true)
%       Check for infiite values?
%
%   --------
%   OUTPUTS:
%   --------
%   - passTest: list(bool) (length changing depending on tests performed)
%       Did 'var_TBC' pass the test.
%       passTest(1) = nan, passTest(2) = 0, passTest(3) = infty, 
%
%   --------
%   ISSUES:
%   --------
%
%   --------
%   TODO:
%   --------
%   - Return a vector of test (?)
%   - Add a masking option based on the convergence status of each pixel.

    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = false;
    end
    
    % Test variable
    check_bool = false;
    
    % Sizes:
    n_layer = length(theta);
    n_scale = zeros(1,n_layer);
    s_image = size(theta{1}.proba{1}(:,:,1));
    n_state = size(theta{1}.proba{1},3);
    
    for layer=1:n_layer
        n_scale(1,layer) = length(theta{layer}.proba);
    end
    
    %% +++ Sanity check:
    for layer=1:n_layer
        
        field = fieldnames(theta{layer});
        for i = 1:numel(field)           
            if not(strcmp(field{i},'distr'))
                tmp_field = theta{layer}.(field{i});
                for scale=1:n_scale(layer)
                    % NaN:
                    nan_c_mat = isnan(tmp_field{scale});
                    % Inf:
                    inf_c_mat = abs(tmp_field{scale}) == inf;

                    if max(max(max(max(nan_c_mat)))) ...
                            || max(max(max(max(inf_c_mat)))) == 1
                        % Update test variable
                        fprintf('+++ check theta: Faulty field: %s at layer %i and scale %i \n',...
                            field{i}, layer, scale)
                        fprintf('+++ check theta: nan_c_mat= %i \n', ...
                            max(max(max(max(nan_c_mat)))))
                        fprintf('+++ check theta: inf_c_mat= %i \n',...
                            max(max(max(max(inf_c_mat)))))

                        check_bool = true;
                        return
                    end
                end
            end
        end
    end
end

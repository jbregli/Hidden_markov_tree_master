function [nan_check_mat, nan_check_bool] = ...
    hmm_Scheck_0nan(var_TBC, var_fnct, var_name, layer, scale, ...
        verbose, nan, zero, infty)
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
    if ~exist('nan','var')
        nan = true;
    end
    if ~exist('zero','var')
        zero = true;
    end
    if ~exist('infty','var')
        infty = true;
    end
    if ~exist('mask','var')
        mask = find(ones(size(var_TBC)));
    end
    
    % Test variable
    nan_check_bool = false * zeros(1,nan + zero + infty);
    
%     % Resize cv_test if necessary:
%     if length(size(cv_test)) ~= length(size(var_TBC))
%         cv_test = repmat(cv_test,1,1,1,size(cv_test,3));
%     end
    
    %% +++ Sanity check:
    if nan
        nan_c_mat = isnan(var_TBC);
        
        if max(max(max(nan_c_mat)))
            % Update test variable
            nan_check_bool(nan) = true;
            
            % Optional print:
            if verbose
                disp([var_fnct ': NAN ' var_name ' at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
        end
    end
    
    if zero
        zero_c_mat = var_TBC == 0;
        
        if max(max(max(zero_c_mat))) == 1
            % Update test variable
            nan_check_bool(nan+zeros) = true;
                        
            % Optional print:
            if verbose
                disp([var_fnct ': 0 ' var_name ' at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
        end
    end
    
    if infty
        inf_c_mat = abs(var_TBC) == inf;
        
        if max(max(max(zero_c_mat))) == 1
            % Update test variable
            nan_check_bool(nan+zeros+infty) = true;
            
            % Optional print:
            if verbose
                disp([var_fnct ': inf ' var_name ' at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
        end
    end
    
    nan_check_bool = max(nan_check_bool);
    nan_check_mat = max(max(nan_c_mat, zero_c_mat), inf_c_mat);
end

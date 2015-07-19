function passTest = hmm_Scheck_0nan(var_TBC, var_fnct, var_name,...
                                            layer, scale, verbose, ...
                                            nan, zero, infty)
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
%   - Add a masking option to ignore some pixels

    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = true;
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
    passTest = true * ones(1,nan + zero + infty);

    %% +++ Sanity check:
    if nan
        if max(max(max(isnan(var_TBC(mask)))))
            % Update test variable
            passTest(nan) = false;
            
            % Optional print:
            if verbose
                disp([var_fnct ': NAN ' var_name ' at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
        end
    end
    if zero
        if max(max(max(var_TBC(mask) == 0))) == 1
            % Update test variable
            passTest(nan+zeros) = false;
            
            % Optional print:
            if verbose
                disp([var_fnct ': 0 ' var_name ' at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
        end
    end
    if infty
        if max(max(max(abs(var_TBC(mask)) == inf))) == 1
            % Update test variable
            passTest(nan+zeros+infty) = false;
            
            % Optional print:
            if verbose
                disp([var_fnct ': inf ' var_name ' at layer ' num2str(layer) ...
                    ' and scale ' num2str(scale)])
            end
        end
    end
    
    % +++ Better if passing a vector but quick fix:
    passTest = min(passTest);
 
end

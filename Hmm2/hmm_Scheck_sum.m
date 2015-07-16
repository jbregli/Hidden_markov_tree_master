function passTest = hmm_Scheck_sum(var_TBC, target, var_fnct, var_name,...
                        target_name, layer, scale, verbose)
% hmm_Scheck_sum: CHECK IF THE GIVEN VARIABLE IS SUMMING TO ITS TARGET.
%   
%   --------
%   INPUTS:
%   --------
%   - var_TBC: Multidimensional Arrays (d=3 or 4)
%       Matrix of size [s_im, n_state] or [s_im, n_state, n_state] holding 
%       the value to check.
%   - target: Multidimensional Arrays (d=2 or 3)
%       Target matrix for the summing of 'ar_TBC' over dimension 3 or 4.
%   - var_fnct: string
%       Function name from where the checking is done. Help for debugging.
%       Display only if 'verbose' is true.
%   - var_name: string
%       'Var_TBC' name. Help for debugging. Display only if 'verbose' is
%       true.
%   - target_name: string
%       'target' name. Help for debugging. Display only if 'verbose' is
%       true.
%   - layer: int
%       Layer number of 'var_TBC'. Help for debugging. Display only if 
%       'verbose' is true.
%   - scale: int
%       Scale number of 'var_TBC'. Help for debugging. Display only if 
%       'verbose' is true.
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%
%   --------
%   OUTPUTS:
%   --------
%   - passTest: bool
%       Did 'var_TBC' pass the test
%       True if the test is validated
%
%   --------
%   ISSUES:
%   --------
%
%   --------
%   IMPROVEMENTS:
%   --------
%
    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = true;
    end

    % Sizes:
    dim = ndims(var_TBC);

    % Test variable
    passTest = true;

    %% +++ Sanity check:
    % dim == 3:
    if dim == 3
        sum_var_TBC = squeeze(sum(var_TBC,3));        
        
        if  any(any(sum_var_TBC - target > 1e-14)) % any(any(sum_var_TBC ~= target))
            % Update test variable
            passTest = false;

            % Optional print:
            if verbose
                disp([var_fnct ': sum_k(' var_name ')' ' != ' target_name ...
                    ' at layer ' num2str(layer) ' scale ' num2str(scale) ])
            end
        end

        % dim == 4:
    elseif dim == 4
        % Number of state:
        n_state = size(var_TBC,3);

        for f_state=1:n_state
            sum_var_TBC = squeeze(sum(var_TBC(:,:,f_state,:),4));

            if any(any(sum_var_TBC - target(:,:,f_state) > 1e-14))  % any(any(sum_var_TBC ~= target(:,:,f_state)))
                % Update test variable
                passTest = false;

                % Optional print:
                if verbose
                    disp([var_fnct ': sum_k(' var_name ') '...
                        '!= ' target_name ' at layer ' num2str(layer) ...
                        ' scale ' num2str(scale) ...
                        ' and for f_state  ' num2str(f_state)])

                    rapport = sum_var_TBC ./ target(:,:,f_state);

                    if all(rapport == rapport(1))
                        disp(['The rapport sum_k(' var_name ') ./ ' ...
                            target_name ' = ' num2str(rapport(1))])
                    end
                end
            end
        end
    else
        disp([var_fnct ' :Dimension of the var_TBC = ' num2str(dim) ...
            ' --> unknown dimension'])
    end
end


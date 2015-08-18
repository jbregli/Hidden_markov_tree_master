function var_TBC = correct_epsilon(var_TBC)
% correct_epsilon: CHECK IF THE GIVEN VARIABLE IS SUMMING TO ITS TARGET.
%   
%   --------
%   INPUTS:
%   --------
%   - var_TBC: Multidimensional Arrays (d=3 or 4)
%       Matrix of size [s_im, n_state] or [s_im, n_state, n_state] holding 
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
%   TODO:
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
    % dim == 2 -- Epsilon uniform
    if dim == 2
        sum_var_TBC = squeeze(sum(var_TBC,2));     
        
        tmp_test = ones(size(sum_var_TBC)) - sum_var_TBC;
        
        % Add the difference to the maximum of the transition values:
        [~, max_ind] = max(var_TBC,[],2);     
        var_TBC(max_ind) = var_TBC(max_ind) + tmp_test;
        
%         if  any(tmp_test ~= 0) % any(any(sum_var_TBC ~= target))
%             % If the difference is small add it to the maximum of the the
%             % transsition values:
%             
%             % Optional print:
%             if verbose
%                 disp([var_fnct ': sum_k(' var_name ')' ' != ' target_name ...
%                     ' at layer ' num2str(layer) ' scale ' num2str(scale) ])
%             end
%         end

    % dim == 4:
    elseif dim == 4
        % Number of state:
        n_state = size(var_TBC,3);

        for f_state=1:n_state
            tmp_var_TBC = squeeze(var_TBC(:,:,f_state,:));
            sum_var_TBC = squeeze(sum(tmp_var_TBC,3));

            tmp_test = ones(size(sum_var_TBC)) - sum_var_TBC;
            
            % Add the difference to the maximum of the transition values:
            [~, max_ind] = max(tmp_var_TBC(:,:,f_state,:),[],3);            
            tmp_var_TBC(max_ind) = tmp_var_TBC(max_ind) + tmp_test;
            
            % Update 'var_TBC':
            var_TBC(:,:,f_state,:) = tmp_var_TBC;

%             
%             if any(any(abs(sum_var_TBC - target(:,:,f_state)) > 1e-15))  % any(any(sum_var_TBC ~= target(:,:,f_state)))
%                 % Update test variable
%                 passTest = false;
% 
%                 % Optional print:
%                 if verbose
%                     disp([var_fnct ': sum_k(' var_name ') '...
%                         '!= ' target_name ' at layer ' num2str(layer) ...
%                         ' scale ' num2str(scale) ...
%                         ' and for f_state  ' num2str(f_state)])
% 
%                     rapport = sum_var_TBC ./ target(:,:,f_state);
% 
%                     if all(rapport == rapport(1))
%                         disp(['The rapport sum_k(' var_name ') ./ ' ...
%                             target_name ' = ' num2str(rapport(1))])
%                     end
%                 end
%             end
        end
        % Update 'var_TBC':
        
        
    else
        disp([var_fnct ' :Dimension of the var_TBC = ' num2str(dim) ...
            ' --> unknown dimension'])
    end
end


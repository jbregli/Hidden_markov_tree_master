function [b_index] = hmm_find_brothers(S, c_layer, c_index)
% hmm_find_parent: FIND THE BROTHERS OF A GIVEN NODE:
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.      
%   - c_layer: int
%       Layer of the child
%   - c_index: int
%       Index of the child
%
%   --------
%   OUTPUTS:
%   --------
%   - b_index: {int}
%       Indexes of the brothers
%
%   --------
%   IMPROVEMENTS:
%   --------

    %% Initialization:
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1, n_layer);
    for l=1:n_layer
        n_elmt(1, l) = length(S{l}.signal);
    end
    
    %% Find parent:
    % Find in the previous layer the index of the father node:
    % Same scale path:
    % The last layer has no children:
    if c_layer == 1
        b_index = [];
    % The first layer has all the second layer as children:
    elseif c_layer == 2
        b_index = (1:n_elmt(1, c_layer));
        % Remove c_index from b_index:
        b_index = b_index(b_index ~= c_index);        
    else
        % Find the father:
        f_layer = c_layer - 1;
        if isfield('parent', S{c_layer}.hmm{c_index})
            f_index = S{c_layer}.hmm{c_index}.parent;
        else
            f_index = hmm_find_parent(S, c_layer, c_index);
        end
        
        % Find all the children of the father
        if isfield('children', S{f_layer}.hmm{f_index})
            b_index = S{f_layer}.hmm{f_index}.parent;
        else
            b_index = hmm_find_children(S, f_layer, f_index);
        end        
        
        % Remove c_index from b_index:
        b_index = b_index(b_index ~= c_index);

        
    end    
        

end


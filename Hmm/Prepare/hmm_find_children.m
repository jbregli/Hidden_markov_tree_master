function [c_index] = hmm_find_children(S, f_layer, f_index)
% hmm_find_parent: FIND THE CHILDREN OF A GIVEN NODE:
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.      
%   - f_layer: int
%       Layer of the father
%   - f_index: int
%       Index of the father
%
%   --------
%   OUTPUTS:
%   --------
%   - c_index: {int}
%       Indexes of the children
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
    if f_layer == n_layer
        c_index = [];
    % The first layer has all the second layer as children:
    elseif f_layer == 1
        c_index = (1:n_elmt(1, f_layer+1));
    else
        % Scale path of the father:
        f_scale =  S{f_layer}.meta.j(1:end-1,f_index);
        % Orientation path of the father:
        f_theta =  S{f_layer}.meta.theta(1:end,f_index);      

        % Find this path and orientation at layer 'l+1':
        tmp_cj = ismember(S{f_layer+1}.meta.j(1:end-2,:)', f_scale', 'rows');
        index_scale = find(tmp_cj == 1);
        % Same orientation
        tmp_ctheta = ismember(S{f_layer+1}.meta.theta(1:end-1,:)', f_theta', 'rows');
        index_theta= find(tmp_ctheta == 1);

        % Find the matching index:
        c_index = [intersect(index_scale, index_theta)];
        
        
    end    
        

end


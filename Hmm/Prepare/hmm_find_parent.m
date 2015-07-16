function [f_index] = hmm_find_parent(S, c_layer, c_index)
% hmm_find_parent: FIND THE PARENT OF A GIVEN NODE:
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
%   - c_index: int
%       Index of the father
%
%   --------
%   IMPROVEMENTS:
%   --------

    %% Initialization:
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1,n_layer);
    for l=1:n_layer
        n_elmt(1,l) = length(S{l}.signal);
    end
    
    %% Find parent:
    % Find in the previous layer the index of the father node:
    % Same scale path:
    if c_layer == 1
        f_index = [];
    elseif c_layer == 1
        f_index = 1;
    else
        % Scale path of the father:
        c_scale = S{c_layer}.meta.j(1:end-2,c_index);
        % Orientation path of the father:
        c_theta = S{c_layer}.meta.theta(1:end-1,c_index);      

        % Find this path and orientation at layer 'l-1':
        tmp_fj = ismember(S{c_layer-1}.meta.j(1:end-1,:)', c_scale', 'rows');
        index_scale = find(tmp_fj == 1);
        % Same orientation
        tmp_ftheta = ismember(S{c_layer-1}.meta.theta', c_theta', 'rows');
        index_theta= find(tmp_ftheta == 1);

        % Find the matching index:
        f_index = intersect(index_scale, index_theta);
    end    
        

end


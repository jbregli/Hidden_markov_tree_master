function [ prun_params] = isfield_Pruning( prun_params)
% isfield_Pruning: CHECK FOR MISSING FIELDS AND ASSIGNED THEM TO DEFAULT
%
%   --------
%   INPUTS:
%   --------
%   - isfield_Pruning: (optional) struct
%        Each fields is a meta parameter for the algo. 
%
%   --------
%   OUTPUTS:
%   --------
%   - prun_params: (optional) struct
%       Each fields is a meta parameter for the pruning.
%           - .rmv_prc: (optional) float(0,1) (default: 0.5)
%               Percentage of highest signal to noise ratio leafs to be
%               removed
%           - .n_iteration: (optional) int (default: 1)
%               Number of iteration over the pruning. If all the children 
%               of a node are removed during the first pruning iteration
%               then the second pass is done on the   ---- TBC


    if ~isfield(prun_params, 'rmv_prc')
        prun_params.n_step = 0.5;
    end
    if ~isfield(prun_params, 'n_iteration')
        prun_params.n_state = 1;
    end    
end


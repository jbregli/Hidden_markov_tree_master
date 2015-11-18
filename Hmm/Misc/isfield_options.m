function [options] = isfield_options( options)
% isfield_EMmeta: CHECK FOR MISSING FIELDS AND ASSIGNED THEM TO DEFAULT
%
%   --------
%   INPUTS:
%   --------
%   - options: (optional) struct
%        Each fields is an option for the algo. 
%
%   --------
%   OUTPUTS:
%   --------
%   - options: (optional) struct
%        Each fields is an option for the algo. 
%           - .verbose: (optional) bool (default= false)
%                   If true then 'hmm_Scheck_sum' displays debugging info.

    if ~isfield(options, 'verbose')
        options.verbose = false;
    end
end


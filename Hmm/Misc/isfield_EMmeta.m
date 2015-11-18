function [ EM_metaparameters] = isfield_EMmeta( EM_metaparameters)
% isfield_EMmeta: CHECK FOR MISSING FIELDS AND ASSIGNED THEM TO DEFAULT
%
%   --------
%   INPUTS:
%   --------
%   - EM_metaparameters: (optional) struct
%        Each fields is a meta parameter for the algo. 
%
%   --------
%   OUTPUTS:
%   --------
%   - EM_metaparameters: (optional) struct
%        Each fields is a meta parameter for the algo. 
%           - .n_step: (optional) int (default: 100)
%               Maximum number of steps for the EM algorithm.
%           - .n_state: (optional) int (default: 2)
%               Number of states for the mixture of models.
%           - .distribution: (optional) str (default: MixtGauss)
%               Distribution to be used for modeling.
%               Options: MixtGauss
%           - .eps_uni: (optional) bool (default: false)
%               If set to true epsilon is considered to be uniform accross 
%               the image (ie: same probability transition from father to
%               child for all the pixel of a same image. 
%           - .mixing: (optional) int (default: 10)
%               Number of step beforetesting convergence.
%           - .cv_sens: (optional) float (default: 1e-5)
%               Minimum increase between two steps to consider convergence.
%           - .cv_steps: (optional) int (default: 7)
%               Number of consecutive steps where the increase has to be 
%               lower than 'cv_sens' to consider convergence.
%           - .cv_ratio: (optional) float[0,1] (default: 0.98) 
%                Ratio of converged pixel before breaking.

    if ~isfield(EM_metaparameters, 'n_step')
        EM_metaparameters.n_step = 100;
    end
    if ~isfield(EM_metaparameters, 'n_state')
        EM_metaparameters.n_state = 2;
    end    
    if ~isfield(EM_metaparameters, 'distribution')
        EM_metaparameters.distribution = 'MixtGauss';
    end    
    if ~isfield(EM_metaparameters, 'eps_uni')
        EM_metaparameters.eps_uni = false;
    end    
    if ~isfield(EM_metaparameters, 'mixing')
        EM_metaparameters.mixing = 10;
    end    
    if ~isfield(EM_metaparameters, 'cv_sens')
        EM_metaparameters.cv_sens = 1e-5;
    end    
    if ~isfield(EM_metaparameters, 'cv_steps')
        EM_metaparameters.cv_steps = 7;
    end    
    if ~isfield(EM_metaparameters, 'cv_ratio')
        EM_metaparameters.cv_ratio = 0.98;
    end    
end


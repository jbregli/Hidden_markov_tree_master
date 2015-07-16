function [alpha, dob] = conditional_DOWN(S, theta, hidStates, cond_up, ...
    verbose)                                                                % OK
% conditional_UP: COMPUTE THE CONDITIONAL DOWN STEP FOR SHMT.
%   
%   See "Statistical inference for Hidden Markov Tree Models and 
%   application to wavelet trees" for more details on the algorithm.
%   See algorithm 4 p12.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%   - hidStates: cell(struct)
%       Cell of structures containing the distributions for S. Obtained 
%       with the function distribution_HS.
%   - cond_up: cell(struct)
%       Cell of structures containing the values of 'Beta' and the 
%       probabilities of the model.
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%
%   --------
%   OUTPUTS:
%   --------
%   - alpha: cell(struct)
%       Cell of structures containing the values of 'Alpha'.
%   - dob: (optional) struct 
%       Debuging Object, place holder to pass along all variable needed for
%       debugging.
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Use sanity checks as stopping signals.

    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = true;
    end

    % Sizes:
    n_layer = length(S);
    n_state = size(theta{1}.proba{1}, 3);
    s_image = size(S{1}.signal{1});

    n_elmt = zeros(1,n_layer);

    % Structure to store the 'Alpha'
    alpha = cell(1, n_layer);

    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);

        % Structure:
        alpha{layer} = cell(1,n_elmt(1,layer));
        % Initialize the matrices:
        for i=1:n_elmt(1,layer)
            alpha{layer}{i} = ones([s_image n_state]);
        end
    end

    %% Initialization:
    % Root of the tree:
    % ---> Root is initialized to 1. This is done in the creation of alpha

    %% Induction:
    for layer=2:n_layer
        for scale=1:n_elmt(1,layer)
            % a_u - alpha
            % a_u(k) = sum_{i=1}^{K}(
            %               a_{rho(u)(i) eps_{u,rho(u)(i,k) B_{rho(u)\u}(i)
            %               P(S_{rho(u)} = i))
            %           / P(S_u = k))

            % Scale and layer of the father node:
            f_layer = layer-1;
            f_scale = S{layer}.hmm{scale}.parent;

            % Find the children 'scale' in the child list:
            [~,loc] = ismember(scale,S{f_layer}.hmm{f_scale}.children);

            alpha{layer}{scale} = ...
                squeeze(sum(...
                repmat(...
                alpha{f_layer}{f_scale} ...                                   % a_{rho(u)(i)
                .* cond_up{f_layer}.beta.excludeChild{f_scale}{loc} ... 	% eps_{u,rho(u)(i,k)
                .* hidStates{f_layer}.ofHiddenStates{f_scale}, ...          % B_{rho(u)\u}(i)
                1,1,1,n_state) ...
                .* theta{layer}.epsilon{scale}, ...                             % P(S_{rho(u)} = i)
                3)) ...
                ./ hidStates{layer}.ofHiddenStates{scale};                      % / P(S_u = k))

            % +++ SANITY CHECKS:
            % nan
            check_nan = hmm_Scheck_0nan(alpha{layer}{scale},...
                'cond_DOWN', 'alpha', layer, scale, verbose);
        end
    end
end


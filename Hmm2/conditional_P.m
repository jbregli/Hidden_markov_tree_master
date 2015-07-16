function [proba, s_check, dob] = conditional_P(S, theta, hidStates, ...
    cond_up, alpha, verbose)
% conditional_P: COMPUTE THE CONDITIONAL PROBABILITIES FOR SHMT.
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
%   - distribution: cell(struct)
%       Cell of structures containing the distributions for S. Obtained
%       with the function distribution_HS.
%   - cond_up: cell(struct)
%       Cell of structures containing the values of 'Beta' and the
%       probabilities of the model.
%   - alpha: cell(struct)
%       Cell of structures containing the values of 'Alpha'.
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%
%   --------
%   OUTPUTS:
%   --------
%   - set_S: cell(struct)
%       Set of structures updated
%   - s_check: (optional) bool
%       Have the conditional probabilities passed the test?
%   - dob: (optional) struct
%       Debuging Object, place holder to pass along all variable needed for
%       debugging.
%
%   --------
%   ISSUES:
%   --------
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Use sanity checks as matricial stopping signals.

    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = true;
    end

    % Sizes and structure to store the distribution:
    n_layer = length(S);
    n_state = size(theta{1}.proba{1}, 3);
    s_image = size(S{1}.signal{1});

    n_elmt = zeros(1,n_layer);

    % Structure to store the 'proba'
    proba = cell(1, n_layer);

    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);

        % Structure:
        proba{layer}.ofNode = cell(1,n_elmt(1,layer));
        proba{layer}.ofNodeAndParent = cell(1,n_elmt(1,layer));

        % Initialize the matrices:
        for i=1:n_elmt(1,layer)
            proba{layer}.ofNode{i} = zeros([s_image n_state]);
            proba{layer}.ofNodeAndParent{i} = zeros([s_image n_state n_state]);
        end
    end

    % Sanity_check:
    s_check = ones(n_layer,max(n_elmt),4);

    %% P(s_u = k | w) & P(s_{rho(u)} = i, s_u = k | w):
    for layer=1:n_layer
        for scale=1:n_elmt(1,layer)
            % P(s_u=k |w) - Proba of hidden state given the node
            % P(s_u=k |w) = a_u(k) B_u(k)
            proba{layer}.ofNode{scale} = alpha{layer}{scale} ...
                .* cond_up{layer}.beta.givenNode{scale};

            if layer > 1
                % P(s_{rho(u)}= i, s_u=k |w) - Proba of hidden state given the node
                % P(s_{rho(u)}= i, s_u=k |w) = B_u(k) eps_{u,rho(u)(i,k)
                %                              a_{rho(u)}(i) B_{rho(u)\u}(i)
                %                              P(s_{rho(u)} = i)
                %                              / P(s_u = k))

                % Scale and layer of the father node:
                f_layer = layer-1;
                f_scale = S{layer}.hmm{scale}.parent;

                % Find the children 'scale' in the child list:
                [~,loc] = ismember(scale,S{f_layer}.hmm{f_scale}.children);

                for f_state=1:n_state                                                   % WARNING: THE FACTOR 2 HERE HAS BEEN ADDED FROM THE CORRECT FORMULA
                    for c_state=1:n_state                                               % IT IS A QUICK PATCH NOT PROVEN YET
                        proba{layer}.ofNodeAndParent{scale}(:,:,f_state,c_state) = ...              % P(s_{rho(u)}= i, s_u=k |w) =
                            cond_up{layer}.beta.givenNode{scale}(:,:,c_state) ...                   % B_{u}(k)
                            ./ hidStates{layer}.ofHiddenStates{scale}(:,:,c_state) ...              % / P(S_u = k))
                            .* theta{layer}.epsilon{scale}(:,:,f_state,c_state) ...                 % eps_{u,rho(u)(i,k)
                            .* alpha{f_layer}{f_scale}(:,:,f_state) ...                             % a_{rho(u)(i)
                            .* cond_up{f_layer}.beta.excludeChild{f_scale}{loc}(:,:,f_state) ...    % B_{rho(u)\u}(i)
                            .* hidStates{f_layer}.ofHiddenStates{f_scale}(:,:,f_state);             % P(S_{rho(u)} = i)
                    end
                end

                %               NOT WORKING
                %                 test_ONP = ...
                %                     repmat(...
                %                     cond_up{layer}.beta.givenNode{scale} ...                    % B_{u}(i)
                %                     ./ hidStates{layer}.ofHiddenStates{scale} ...               % / P(S_u = k))
                %                     .* cond_up{f_layer}.beta.excludeChild{f_scale}{loc} ...     % B_{rho(u)\u}(i)
                %                     .* alpha{f_layer}{f_scale} ...                              % a_{rho(u)(i)
                %                     .* hidStates{f_layer}.ofHiddenStates{f_scale}, ...          % P(S_{rho(u)} = i)
                %                     1,1,1,n_state) ...
                %                     .* theta{layer}.epsilon{scale};                                 % eps_{u,rho(u)(i,k)
                %
                %                 if test_ONP ~= proba{layer}.ofNodeAndParent{scale}
                %                     disp('fail')
                %                 else
                %                     disp('yattaaa?')
                %                 end

                % +++ SANITY CHECKS:
                % 0 or Nan:
                s_check(layer,scale,1) = ...
                    hmm_Scheck_0nan(proba{layer}.ofNodeAndParent{scale},...
                    'cond_P', 'proba_of_node_and_parents',...
                    layer, scale, verbose);
                % Sum:
                s_check(layer,scale,2) = ...
                    hmm_Scheck_sum(proba{layer}.ofNodeAndParent{scale},...
                    proba{f_layer}.ofNode{f_scale}, ...
                    'Cond_P', 'P(s_{rho(u)}= i, s_u=k |w)',...
                    'P(s_{rho(u)}= i)', layer, scale, verbose);
            end
            % +++
            s_check(layer,scale,3) = ...
                hmm_Scheck_0nan(proba{layer}.ofNode{scale},...
                'cond_P', 'proba_of_node',...
                layer, scale, verbose);

            % sum:
            s_check(layer,scale,4) = ...
                hmm_Scheck_sum(proba{layer}.ofNode{scale}, ...
                ones(size(sum(proba{layer}.ofNode{scale},3))),...
                'Cond_P', ' P(s_u=k |w)', '[1]', layer, scale, verbose);
            % +++
        end
    end

    s_check = min(min(min(s_check)));
end


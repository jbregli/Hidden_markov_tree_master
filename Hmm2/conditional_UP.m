function [cond_up, dob] = conditional_UP(S, theta, hidStates, verbose)      % OK
% conditional_UP: COMPUTE THE CONDITIONAL UP STEP FOR SHMT.
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
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%
%   --------
%   OUTPUTS:
%   --------
%   - cond_up: cell(struct)
%       Cell of structures containing the values of 'Beta' and the 
%       probabilities of the model.
%   - dob: (optional) struct 
%       Debuging OBject, place holder to pass along all variable needed for
%       debugging.
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Add threshold for P_{theta_k}(w_u)} in the arguments (optional maybe)
%   - Use sanity checks as stopping signals.
%   - Check for unused fields in the structure 'cond_up' --> how is l used
%     in the next steps?

    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = true;
    end

    % Threshold for 'model_proba':
    mprob_thres = 1e-3;
    min_thres = NaN; %1e-3;
    max_thres = NaN; %100;

    % Sizes:
    n_layer = length(S);
    n_state = size(theta{1}.proba{1}, 3);
    s_image = size(S{1}.signal{1});

    n_elmt = zeros(1,n_layer);

    % Structure to store the 'Beta'
    cond_up = cell(1, n_layer);

    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);

        % Structure:
        cond_up{layer}.Pvalue = cell(1,n_elmt(1,layer));
        cond_up{layer}.beta.givenNode = cell(1,n_elmt(1,layer));
        cond_up{layer}.beta.givenParent = cell(1,n_elmt(1,layer));
        cond_up{layer}.beta.excludeChild = cell(1,n_elmt(1,layer));
        % Initialization of the matrices:
        for i=1:n_elmt(1,layer)
            cond_up{layer}.Pvalue{i} = zeros([s_image n_state]);
            cond_up{layer}.beta.givenNode{i} = zeros([s_image n_state]);
            cond_up{layer}.beta.givenParent{i} = zeros([s_image n_state]);
            cond_up{layer}.beta.excludeChild{i} = cell(1,length(S{layer}.hmm{i}.children));
            for j=1:length(S{layer}.hmm{i}.children)
                cond_up{layer}.beta.excludeChild{i}{j} = zeros([s_image n_state]);
            end
        end
    end

    %% Initialisation:
    % P_{theta_k}(w_u)}
    for layer=n_layer:-1:1
        for scale=1:n_elmt(1,layer)
            cond_up{layer}.Pvalue{scale} = model_proba(S, theta, layer, ...
                scale, mprob_thres);

            % +++
            check_nan = hmm_Scheck_0nan(cond_up{layer}.Pvalue{scale},...
                'cond_UP_init', 'beta_Pvalue',...
                layer, scale, verbose);
            % +++
        end
    end

    % Loop over the leaves of the tree:
    for layer=n_layer:-1:1
        for scale=1:n_elmt(1,layer)
            % A leave has no children:
            if isempty(S{layer}.hmm{scale}.children)
                %---------------------------------------------------------%
                % B_u - Beta givenNode:                     LEAFS         %
                % B_u(k) = P_{theta_k}(w_u)} P(S_u=k)                     %
                %            / sum_{i=1}^{K} (P_{theta_i}(w_u)} P(S_u=i)) %
                %---------------------------------------------------------%
                % a) Compute the numerator
                % b) and the sum for the denominator a:
                cond_up{layer}.beta.givenNode{scale} = ...
                    cond_up{layer}.Pvalue{scale} ...
                    .* hidStates{layer}.ofHiddenStates{scale};

                tmp_sumBeta = sum(cond_up{layer}.beta.givenNode{scale},3);

                cond_up{layer}.beta.givenNode{scale} = ...
                    cond_up{layer}.beta.givenNode{scale} ...
                    ./ repmat(tmp_sumBeta,1,1,n_state);

                % +++ SANITY CHECKS:
                % Beta given Node should sum to [1]:
                check_sum = hmm_Scheck_sum(cond_up{layer}.beta.givenNode{scale}, ...
                    ones(size(sum(cond_up{layer}.beta.givenNode{scale},3))),...
                    'Cond_UP', 'beta_givenNode', '[1]', layer, scale, verbose);
                % +++

                %---------------------------------------------------------%
                % B_{u,rho(u)} - Beta givenParent:          LEAFS         %
                % B_{u,rho(u)}(j) = sum_{i=1}^{K}(                        %
                %                               B_u(i) eps(j,i)/P(s_u=i)) %
                %                           .P(s_{rho(u)}=j)              %
                %---------------------------------------------------------%
                % Scale and layer of the father node:
                f_layer = layer-1;
                f_scale = S{layer}.hmm{scale}.parent;

                for f_state=1:n_state                                       % ----> OK
                    cond_up{layer}.beta.givenParent{scale}(:,:,f_state) = ...
                        sum(cond_up{layer}.beta.givenNode{scale} ...
                        .* squeeze(theta{layer}.epsilon{scale}(:,:,f_state,:)) ...
                        ./ hidStates{layer}.ofHiddenStates{scale}, ...
                        3) ...
                        .* hidStates{f_layer}.ofHiddenStates{f_scale}(:,:,f_state);
                end

                % +++ SANITY CHECKS:
                % Beta given Parents should sum to [1]:
                check_sum = hmm_Scheck_sum(cond_up{layer}.beta.givenParent{scale}, ...
                    ones(size(sum(cond_up{layer}.beta.givenParent{scale},3))),...
                    'Cond_UP', 'beta_givenParent', '[1]', layer, scale, verbose);
                % +++

                % Optional thresholding:
                if not(isnan(min_thres))
                    cond_up{layer}.beta.givenParent{scale} = ...
                        cond_up{layer}.beta.givenParent{scale} ...
                        .* (cond_up{layer}.beta.givenParent{scale}>min_thres) ...
                        + min_thres ...
                        .* (cond_up{layer}.beta.givenParent{scale}<=min_thres);
                end

                %---------------------------------------------------------%
                % l_u:                                      LEAFS         %
                %---------------------------------------------------------%
                cond_up{layer}.l{scale} = 0;

                % +++
                check_nan = hmm_Scheck_0nan(cond_up{layer}.beta.givenNode{scale},...
                    'cond_UP_init', 'beta_given_Node',...
                    layer, scale, verbose);
                check_sum = hmm_Scheck_0nan(cond_up{layer}.beta.givenParent{scale},...
                    'cond_UP_init', 'beta_given_parents',...
                    layer, scale, verbose);
            end
        end
    end

    %% Induction:
    % Bottom-Up loop on the nodes of the tree
    for layer=(n_layer-1):-1:1
        for scale=1:n_elmt(1,layer)
            % Make sure this node is not a leaf:
            if not(isempty(S{layer}.hmm{scale}.children))

                n_children = length(S{layer}.hmm{scale}.children);

                %---------------------------------------------------------%
                % M_u:                                      NODES         %
                % M_u = sum_{i=1}^{K} (                                   %
                %            P_{theta_i}(w_u)} Prod_{t \in c(u)}(B_{t,u}(i))
                %              / P(S_u=i) ^{n_u-1}                        %
                %---------------------------------------------------------%
                tmp_prodBeta = ones(size(cond_up{layer}.beta.givenParent{scale}));
                tmp_sumL = 0;
                % Loop over the children of the node:
                for t=1:n_children
                    % Scale and layer of the child node:
                    c_layer = layer + 1;
                    c_scale = S{layer}.hmm{scale}.children(t);

                    % Product over the children for M_u and B_u
                    tmp_prodBeta = tmp_prodBeta .* ...
                        cond_up{c_layer}.beta.givenParent{c_scale};

                    % Sum over the children for l_u
                    tmp_sumL = tmp_sumL + cond_up{c_layer}.l{c_scale};
                end

                M_u = sum(cond_up{layer}.Pvalue{scale} .* tmp_prodBeta ...
                    ./ (hidStates{layer}.ofHiddenStates{scale} ...
                    .^ (n_children-1)),...
                    3);

                %---------------------------------------------------------%
                % l_u:                                      NODES         %
                % l_u = log(M_u) + sum_{t \in c(u)}(l_{t})                %
                %---------------------------------------------------------%
                cond_up{layer}.l{scale} = log(M_u) + tmp_sumL;

                %---------------------------------------------------------%
                % B_u - Beta givenNode:                     NODES         %
                % B_u(k) = P_{theta_k}(w_u)} Prod_{t \in c(u)}(B_{t,u}(k))%
                %              / (P(S_u=k) ^{n_u-1} M_u)                  %
                %---------------------------------------------------------%
                cond_up{layer}.beta.givenNode{scale} = ...
                    cond_up{layer}.Pvalue{scale} .* tmp_prodBeta ...
                    ./ (hidStates{layer}.ofHiddenStates{scale} ...
                    .^ (n_children-1) ...
                    .* repmat(M_u,1,1,n_state));

                % +++ SANITY CHECKS:
                % sum:
                check_sum = hmm_Scheck_sum(cond_up{layer}.beta.givenNode{scale}, ...
                    ones(size(sum(cond_up{layer}.beta.givenNode{scale},3))),...
                    'Cond_UP', 'beta_givenNode', '[1]', layer, scale, verbose);
                % +++

                %----------------------------------------------------------%
                % B_{u\c(u)} - Beta exclude child:           NODES         %
                % B_{u\c(u)}(k) = B_u(k) / B_{t,u}(k)                      %
                %----------------------------------------------------------%
                % Layer of the child nodes:
                c_layer = layer + 1;
                for t=1:n_children
                    % Scale of this child node:
                    c_scale = S{layer}.hmm{scale}.children(t);

                    cond_up{layer}.beta.excludeChild{scale}{t} = ...
                        cond_up{layer}.beta.givenNode{scale} ...
                        ./ cond_up{c_layer}.beta.givenParent{c_scale};

                    % Optional thresholding:
                    if not(isnan(max_thres))
                        cond_up{layer}.beta.excludeChild{scale}{t} = ...
                            cond_up{layer}.beta.excludeChild{scale}{t} ...
                            .* (cond_up{layer}.beta.excludeChild{scale}{t}<=max_thres) ...
                            + max_thres ...
                            .* (cond_up{layer}.beta.excludeChild{scale}{t}>max_thres);
                    end

                    %                   if max(max(max(cond_up{layer}.beta.excludeChild{scale}{t} >100)))==1
                    %                       disp(['cond_UP: beta_exculde_child too big at ' ...
                    %                             ' layer ' num2str(layer) ' and scale ' ...
                    %                             num2str(scale)])
                    %                   end
                end

                %---------------------------------------------------------%
                % B_{u,rho(u)} - Beta givenParent:          NODES         %
                % B_{u,rho(u)} = sum_{i=1}^{K} (B_u(i) eps_{u,rho(u)}(j,i)%
                %                                    / P(s_u=i))          %
                %                           .P(s_{rho(u)}=j)              %
                %---------------------------------------------------------%
                if layer > 1
                    % Scale and layer of the father node:
                    f_layer = layer-1;
                    f_scale = S{layer}.hmm{scale}.parent;
                    for f_state=1:n_state
                        cond_up{layer}.beta.givenParent{scale}(:,:,f_state) = ...
                            sum(cond_up{layer}.beta.givenNode{scale} ...
                            .* squeeze(theta{layer}.epsilon{scale}(:,:,f_state,:)) ...
                            ./ hidStates{layer}.ofHiddenStates{scale}, ...
                            3) ...
                            .* hidStates{f_layer}.ofHiddenStates{f_scale}(:,:,f_state);
                    end

                    % +++ SANITY CHECKS:
                    % sum:
                    check_sum = hmm_Scheck_sum(cond_up{layer}.beta.givenParent{scale}, ...
                        ones(size(sum(cond_up{layer}.beta.givenParent{scale},3))),...
                        'Cond_UP', 'beta_givenParent', '[1]', layer, scale);
                    % +++

                    %                   % +++
                    %                   if min(min(min(cond_up{layer}.beta.givenParent{scale} <=0.001)))==1
                    %                       disp(['cond_UP: beta_given_parents too small at layer ' num2str(layer) ' and scale ' num2str(scale)])
                    %                   end
                end
            end
        end
    end
end


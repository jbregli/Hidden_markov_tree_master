function [cond_up, dob] = hmm_MAP(S_inp, theta, verbose)      % OK
% hmm_MAP: COMPUTE THE MAP OF AN IMAGE.
%   
%   See "Statistical inference for Hidden Markov Tree Models and 
%   application to wavelet trees" for more details on the algorithm.
%   See algorithm 4 p11 & p14.
%
%   The aim of the MAP algorithm is to find the opti;al hidden tree
%   \hat{h_1}=(\hat{s_1}...\hat{s_n}) maximizing P(H_1=h_1|T_1) and the
%   value \hat{P} of t he maximum.
%
%   --------
%   INPUTS:
%   --------
%   - S_input: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.
%       Input to be classified.
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib used
%       to store the model parameters.
%
%   --------
%   OUTPUTS:
%   --------
%   - P_hat: float
%       probability of the sequence given the observed data and the model.
%   - H_tree: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib used
%       to store the optimal states sequence.  
%   - dob: (optional) struct 
%       Debuging OBject, place holder to pass along all variable needed for
%       debugging.
%
%   --------
%   TODO:
%   --------
%   - All

    %% Preparation:
    % Arguments:
    if ~exist('verbose','var')
        verbose = false;
    end

    % Sizes:
    n_layer = length(S_inp);
    n_state = size(theta{1}.proba{1}, 3);
    s_image = size(S_inp{1}.signal{1});

    n_scale = zeros(1,n_layer);

    % Structure to store 'H_tree' and Pvalue (tmp):
    H_tree = cell(1, n_layer);
    variable = cell(1, n_layer);

    for layer=1:n_layer
        n_scale(1,layer) = length(S_inp{layer}.signal);

        % Structure:
        H_tree{layer} = cell(1,n_scale(1,layer));
        variable{layer}.Pvalue = cell(1,n_scale(1,layer));
        variable{layer}.delta = cell(1,n_scale(1,layer));
        variable{layer}.gamma = cell(1,n_scale(1,layer));

    end

    %% Preliminary:
        % P_{theta_k}(w_u)}
    for layer=n_layer:-1:1
        for scale=1:n_scale(1,layer)
            variable{layer}.Pvalue{scale} = ...
                model_proba(S_inp, theta, layer, scale);

            % +++
            check_nan = hmm_Scheck_0nan(variable{layer}.Pvalue{scale},...
                'Hmm_MAP', 'P(w_u)', layer, scale, verbose);
            % +++
            
            % Delta:
            % Node
            % NodeAndParent
            
            % Gamma
            
        end
    end
    
    %% Initialisation:
    % Loop over the leaves of the tree:
    for layer=n_layer:-1:1
        for scale=1:n_scale(1,layer)
            % A leave has no children:
            if isempty(S_inp{layer}.hmm{scale}.children)
                %---------------------------------------------------------%
                % D_u - Delta :                             LEAFS         %
                % D_u(k) = B_u(k)                                         %
                %        = P_{theta_k}(w_u)} P(S_u=k)                     %
                %            / sum_{i=1}^{K} (P_{theta_i}(w_u)} P(S_u=i)) %
                %---------------------------------------------------------%
                % a) Compute the numerator
                % b) and the sum for the denominator a:
                
                %%%%%%%%%%%%%%%% START FROM HERE (NEED TO DEFINE SOME
                %%%%%%%%%%%%%%%% VARIABLES)
                
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
                f_scale = S_inp{layer}.hmm{scale}.parent;

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
        for scale=1:n_scale(1,layer)
            % Make sure this node is not a leaf:
            if not(isempty(S_inp{layer}.hmm{scale}.children))

                n_children = length(S_inp{layer}.hmm{scale}.children);

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
                    c_scale = S_inp{layer}.hmm{scale}.children(t);

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
                    c_scale = S_inp{layer}.hmm{scale}.children(t);

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
                    f_scale = S_inp{layer}.hmm{scale}.parent;
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


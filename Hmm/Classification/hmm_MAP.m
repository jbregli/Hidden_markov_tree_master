function [P_hat, H_tree, dob] = hmm_MAP(S_inp, theta, verbose)
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
    var = cell(1, n_layer);

    for layer=1:n_layer
        n_scale(1,layer) = length(S_inp{layer}.signal);

        % Structure:
        H_tree{layer} = cell(1,n_scale(1,layer));
        var{layer}.Pvalue = cell(1,n_scale(1,layer));
        var{layer}.delta.node = cell(1,n_scale(1,layer));
        var{layer}.delta.nodeAndParents = cell(1,n_scale(1,layer));
        var{layer}.gamma = cell(1,n_scale(1,layer));

    end

    %% Initialisation:
    % Hidden state probabibility:
    hidStates = conditional_HIDDEN(S_inp, theta, verbose);
    % Up phase:
    cond_up= conditional_UP(S_inp, theta, hidStates, verbose);
    
    for layer=n_layer:-1:1
        for scale=1:n_scale(1,layer)
            % A leave has no children:
            if isempty(S_inp{layer}.hmm{scale}.children)            
                %-------------------------------------------------------------%
                % D_u(k):                                        LEAF         %
                % D_u(k) = B_u(k)                                             %
                %-------------------------------------------------------------%
                var{layer}.delta.node{scale} = ...
                    cond_up{layer}.beta.givenNode{scale};

                %-------------------------------------------------------------%
                % D_{u,\rho(u)}(j):                              LEAF         %
                % D_{u,\rho(u)}(j) =                                          %
                %       max_{i \in states}(D_u(i) eps(j,i) / P(S_u=i))        %
                %           * P(S_{\rho(u)}=j)                                %
                % G_u(j):                                        LEAF`        %
                % G_u(j) = argmax_{i \in states}(D_u(i)eps(j,i)/P(S_u=i))     %
                %-------------------------------------------------------------%
                % Initialize gamma at 1 and Delta_nodeAndParents at 0:
                var{layer}.gamma{scale} = ones(s_image);
                var{layer}.delta.nodeAndParents{scale} = zeros([s_image n_state]);
                
                % Maximum over states:
                for f_state=1:n_state
                    % Scale and layer of the father node:
                    f_layer = layer-1;
                    f_scale = S_inp{layer}.hmm{scale}.parent;

                    % Maximum over child's states:
                    tmp_prod =  var{layer}.delta.node{scale} ...
                        .* squeeze(theta{layer}.epsilon{scale}(:,:,f_state,:)) ...
                        ./ hidStates{layer}.ofHiddenStates{scale};

                    tmp_max = tmp_prod(:,:,1);
                    for c_state=2:n_state
                        tmp_max = bsxfun(@max,tmp_max,tmp_prod(:,:,c_state));

                        % Update argmax:
                        var{layer}.gamma{scale}(...
                            tmp_max == tmp_prod(:,:,c_state)) = c_state;
                    end

                    % D_{u,\rho(u)}(j):
                     var{layer}.delta.nodeAndParents{scale}(:,:,f_state) = ...
                        tmp_max ...
                        .* hidStates{f_layer}.ofHiddenStates{f_scale}(:,:,f_state);    
                end
            end
        end
    end
    
    %% Induction:
    % Bottom-Up loop on the nodes of the tree
    for layer=(n_layer-1):-1:1
        for scale=1:n_scale(1,layer)
            % Make sure this node is not a leaf:
            if not(isempty(S_inp{layer}.hmm{scale}.children))
                % Number of children:
                n_children = length(S_inp{layer}.hmm{scale}.children);
                
                %---------------------------------------------------------%
                % D_u(k):                                      NODES      %
                % D_u(k) = B_u(k)                                         %
                %        = P_{theta_k}(w_u) *                             %
                %           Prod_{t \in c(u)}(D_{t,u}(k))                 %
                %              / (M_u * P(S_u=k)^{n_u-1})                 %
                %---------------------------------------------------------%
                tmp_prodDelta = ...
                    ones([s_image n_state]);

                % Loop over the children of the node:
                for t=1:n_children
                    % Scale and layer of the child node:
                    c_layer = layer + 1;
                    c_scale = S_inp{layer}.hmm{scale}.children(t);

                    % Product over the children for D_u:
                    tmp_prodDelta = tmp_prodDelta .* ...
                        var{c_layer}.delta.nodeAndParents{c_scale};
                end

                var{layer}.delta.node{scale} = ...
                        cond_up{layer}.Pvalue{scale} .* tmp_prodDelta ...
                        ./ (hidStates{layer}.ofHiddenStates{scale} ...
                        .^ (n_children-1) ...
                        .* repmat(cond_up{layer}.M{scale},1,1,n_state));

                % +++ SANITY CHECKS:
                % sum:
                check_sum = ...
                    hmm_Scheck_sum(var{layer}.delta.node{scale}, ...
                    ones(size(sum(var{layer}.delta.node{scale},3))),...
                    'MAP', 'delta_givenNode', '[1]', layer, scale, verbose);
                % +++

                %---------------------------------------------------------%
                % D_{u,\rho(u)}(j):                            NODES      %
                % D_{u,\rho(u)}(j) =                                      %
                %       max_{i \in states}(D_u(i) eps(j,i) / P(S_u=i))    %
                %           * P(S_{\rho(u)}=j)                            %
                % G_u(j):                                      NODES      %
                % G_u(j) = argmax_{i \in states}(D_u(i)eps(j,i)/P(S_u=i)) %
                %---------------------------------------------------------%
                % Initialize gamma at 1 and Delta_nodeAndParents at 0:
                var{layer}.gamma{scale} = ones(s_image);
                var{layer}.delta.nodeAndParents{scale} = zeros([s_image n_state]);
                
                % Maximum over states:
                for f_state=1:n_state
                    % Scale and layer of the father node:
                    f_layer = layer-1;
                    f_scale = S_inp{layer}.hmm{scale}.parent;
                    
                    % Maximum over child's states:
                    tmp_prod =  var{layer}.delta.node{scale} ...
                        .* squeeze(theta{layer}.epsilon{scale}(:,:,f_state,:)) ...
                        ./ hidStates{layer}.ofHiddenStates{scale};
                    
                    tmp_max = tmp_prod(:,:,1);
                    for c_state=2:n_state
                        tmp_max = bsxfun(@max,tmp_max,tmp_prod(:,:,c_state));
                        
                        % Update argmax:
                        var{layer}.gamma{scale}(...
                            tmp_max == tmp_prod(:,:,c_state)) = c_state;
                    end
                    
                    % D_{u,\rho(u)}(j):
                    if layer > 1
                         var{layer}.delta.nodeAndParents{scale}(:,:,f_state) = ...
                            tmp_max ...
                            .* hidStates{f_layer}.ofHiddenStates{f_scale}(:,:,f_state);    
                    end
                end
                
                % +++ SANITY CHECKS:
                % sum:
                check_sum = hmm_Scheck_sum(var{layer}.delta.nodeAndParents{scale}, ...
                    ones(size(sum(var{layer}.delta.nodeAndParents{scale},3))),...
                    'MAP', 'delta_givenParent', '[1]', layer, scale, verbose);
                % +++
            end
        end
    end
    
    %% Termination:
    % \hat(P) = max(P(H_1|T_1))
    % \hat(P) = max_{i \in states}(D_u(i))
    P_hat = var{1}.delta.node{1}(:,:,1);
    H_tree{1}{1} = ones(s_image);
    
    for state=2:n_state
        P_hat = bsxfun(@max, P_hat, var{1}.delta.node{1}(:,:,state));
        
        % Update argmax:
        H_tree{1}{1}(P_hat == var{1}.delta.node{1}(:,:,state)) = state;
    end   
    
    %% Downward tracking:
    % Creation of the tree from the root
    for layer=2:n_layer
        for scale=1:n_scale(1,layer)
            % Scale and layer of the father node:
            f_layer = layer-1;
            f_scale = S_inp{layer}.hmm{scale}.parent;
            
            % Optimal state:
            H_tree{layer}{scale} = ...
                var{layer}.gamma{scale}(H_tree{f_layer}{f_scale});
        end
    end
end


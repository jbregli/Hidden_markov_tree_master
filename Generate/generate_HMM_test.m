function [hmm_net, names] = generate_HMM(S, n_state)
% generate_HMM: Create a bayesian network to descr

% Hidden nodes: 
% Discrete node encoding the state of the ST coef. 
% They are connected to:
%   - The observed continuous node encoding the ST coef distribution
%   - The unobserved discrete nodes encoding the state of the ST coef of the
%     children of this node (in the ST tree)
%
% Visible node:
% Continuous node encoding the ST coef distribution.
% They are connected only to their state.
    
    %% Initialization:
    % Arguments:
    if nargin < 2
        n_state = 2;
    end
    
    % Sizes:
    n_layer = length(S);   
    
    % Number of nodes required:
    n_visible = 0;
    
    % Visible nodes - ST coef
    n_scale = cell(n_layer,1);
    for l=1:n_layer
        n_scale{l} = length(S{l}.signal);
        n_visible = n_visible + length(S{l}.signal);
    end
    % Hidden nodes - States:
    n_hidden = n_visible;
    
    N = n_visible + n_hidden;
    %disp(['+++ N = ' num2str(N)])

    % Names
    names = cell(1,N);
    j = 1;
    for f_layer=1:(n_layer)
        % Loop over all the node in this layer:
        for f_index=1:n_scale{f_layer}    
            names{j} = ['S(' num2str(f_layer) '-' num2str(f_index) ')'];
            names{n_hidden + j} = ['s(' num2str(f_layer) '-' num2str(f_index) ')'];
            j = j + 1;
        end
    end
    
    %% Direct Acyclique Graph:
    % Indexing:
    %   Hidden node: 1:n_hidden
    %   Visible node: (n_hidden+1):N
    dag = zeros(N,N);
    
    % Connections H-V --> a hidden node H_i of index i is connected to the
    % visible node V_i of index (n_hidden+i).
    for i=1:1:n_hidden
        dag(i, n_hidden + i) = 1;
    end
    
    % Connections H-H --> a hidden node H_i of index i is connected to its
    % children hidden node
    c_offset = 1;
    
    % Loop over the layers:
    % (exept the last one which has no visible child)
    for f_layer=1:(n_layer-1)
        % Loop over all the node in this layer:
        for f_index=1:n_scale{f_layer}
            % Find path to the father:
            f_scale = S{f_layer}.meta.j(:,f_index)';
            f_theta = S{f_layer}.meta.theta(:,f_index);

            % Find children:
            f_l_scale = length(S{f_layer}.meta.j(:,f_index));
            f_l_theta = length(S{f_layer}.meta.theta(:,f_index));
            c_l_scale = length(S{f_layer+1}.meta.j(:,1));    
            c_l_theta = length(S{f_layer+1}.meta.theta(:,1));

            % Find in the next layer the index of the children of the father node:
            % Same scale path:
            if f_layer == 
                % all the nodes of the 2nd layers are children:
                c_index = 1:length(S{f_layer+1}.meta.j);
            else
                temp = ismember(S{f_layer+1}.meta.j((c_l_scale-(f_l_scale-1)):c_l_scale,:)',...
                                                                  f_scale, 'rows');
                index_scale = find(temp == 1);
                % Same orientation
                temp = ismember(S{f_layer+1}.meta.theta((c_l_theta-(f_l_theta-1)):c_l_theta,:)',...
                                                                  f_theta', 'rows');
                index_theta= find(temp == 1);
                % Find the matching index:
                c_index = intersect(index_scale, index_theta);
            end
            
            % Update the connection graph:
            dag(f_index,c_offset + c_index)=1;
        end
        c_offset = c_offset + length(c_index);
    end
        
    % +++ Test
    % draw_graph(dag);
    
    %% Bayes net shell:
    % In addition to specifying the graph structure, we must specify the size
    % and type of each node.
    %  If a node is DISCRETE, its size is the number of possible values each
    % node can take on.
    % if a node is CONTINUOUS, it can be a vector, and its size is the length
    % of this vector.
    discrete_nodes = 1:n_hidden;
    node_sizes = n_state .* ones(1,n_hidden);
    
    hmm_net = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes,...
                                                           'names', names);
end


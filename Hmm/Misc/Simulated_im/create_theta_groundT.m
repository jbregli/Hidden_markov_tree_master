function [ theta_gt ] = create_theta_groundT(S, dist_param, epsilon, ...
    rn_prob)
% create_theta_groundT:  CREATE THE STRUCTURE TO STORE THE TREE TRUE 
%   PARAMETERS.
%   Given a scattering transform (cell of structure), this functions 
%   associates to each node the required fields for the HMM modeling.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Set of structures obtained with the function 'scat' of the 
%       'scatnet' lib.      
%   - dist_param: cell(float)
%       Matrix giving the mean and variance of each gaussian to be mixed.
%   - epsilon: matrix(float)
%       Matrix giving the transition probabilities
%   - rn_prob: matrix(float)
%       Matrix giving the initial probability for the state of the root
%       node
%
%   --------
%   OUTPUTS:
%   --------
%   - theta_gt: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib.
%       theta{layer}.proba{index}
%       theta{layer}.epsilon{index}
%       theta{layer}.mu{index}
%       theta{layer}.sigma{index}
%       theta{layer}.distr{index}
%
%   --------
%   IMPROVEMENTS:
%   --------
%

    %% Initialization:   
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);
    end
    s_im = size(S{1}.signal{1});
    n_state = length(rn_prob);
    
    % Theta:
    theta_gt = cell(1,n_layer);
       
    %% Theta field:
    % Go through S to create the correct data structure:
    % Over all the layers
    for layer=1:n_layer
        % Create the hmm field - a cell of same size as S.signal:
        theta_gt{layer}.proba = cell(1, n_elmt(1,layer));
        theta_gt{layer}.epsilon = cell(1, n_elmt(1,layer));
        theta_gt{layer}.mu = cell(1, n_elmt(1,layer));
        theta_gt{layer}.sigma = cell(1, n_elmt(1,layer));
        theta_gt{layer}.distr = 'MixtGauss';
        
        % Over all the orientation and scale:
        for scale=1:n_elmt(1,layer)
            
            % PROBABILITY:
            % First layer: Multinomial distribution - Tabular CPD :
            if layer ==1
                % P_{S1}(m):
                theta_gt{layer}.proba{scale} = ones([s_im, n_state]);
                for m=1:n_state
                    theta_gt{layer}.proba{scale}(:,:,m) = rn_prob(m) .* theta_gt{layer}.proba{scale}(:,:,m);
                end
            else
                theta_gt{layer}.proba{scale} = zeros([s_im, n_state]);
            end
            
            % MEAN MUS
            % Mean parameter - as many mean params as states
            theta_gt{layer}.mu{scale} = ones([s_im, n_state]);
            for m=1:n_state
                theta_gt{layer}.mu{scale}(:,:,m) = dist_param{m}{1} .* theta_gt{layer}.mu{scale}(:,:,m);
            end           
            
            % VARIANCE SIGMA parameter - as many var params as states
            theta_gt{layer}.sigma{scale} = ones([s_im, n_state]);
            for m=1:n_state
                theta_gt{layer}.sigma{scale}(:,:,m) = dist_param{m}{2} .* theta_gt{layer}.sigma{scale}(:,:,m);
            end
            
            
            % TRANSITION PROBABILITIES - Tabular CPD
            % For each pixel n_state x n_state posible transitions stored
            % in the 3rd dimension of the matrix using matlab single
            % indexing convention.
            % The first state (rows) describes the father node, the second
            % state (colums) describe the "studied" child node.
            % e.g.: n_state = 2
            % | (1)L-L  (3)L-H | ---> [(1-11)L-L  (2-21)H-L  (3-12)L-H  (4-22)H-H]
            % | (2)H-L  (4)H-H |
            % To access to L-H at layer l and scale j:
            % S{l}.hmm{j}.hidden.epsilon(:,:,sub2ind([n_state,n_state], 1, 2))
            theta_gt{layer}.epsilon{scale} = ones([s_im, n_state , n_state]);
            for m=1:n_state
                for n=1:n_state
                    theta_gt{layer}.epsilon{scale}(:,:,m,n) = epsilon(m,n) .* theta_gt{layer}.epsilon{scale}(:,:,m,n);
            end
        end
    end
end



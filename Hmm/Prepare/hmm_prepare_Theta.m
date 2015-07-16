function [theta] = hmm_prepare_Theta(set_S, n_state, distribution, ...
                                     eps_uni, verbose)
% hmm_prepare_Theta:  CREATE THE STRUCTURE TO STORE THE TREE PARAMETERS.
%   Given a scattering transform (cell of structure), this functions 
%   associates to each node the required fields for the HMM modeling.
%
%   --------
%   INPUTS:
%   --------
%   - set_S: cell(cell(struct))
%       Set of structures obtained with the function 'scat' of the 
%       'scatnet' lib.      
%   - n_state: (optional) int (default: 2)
%       Number of states the values of the scattering coefficient will be
%       partitioned into.
%       eg: n_state = 2 --> Small (1) and High (2)
%           n_state = 3 --> Small (1), Medium (2) and High (3)...
%   - distribution: (optional) str (default: MixtGauss)
%       Distribution to be used for modeling
%       Options: MixtGauss
%   - eps_uni: (optional) bool (default: true)
%       If set to true epsilon is considered to be uniform accross the
%       image (ie: same probability transition from father to child for all
%       the pixel of a same image
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%
%   --------
%   OUTPUTS:
%   --------
%   - theta: cell(struct)
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
%   - possibility to include priors
%   - use a vargin method
%   - add assert for control
%   - Meta parameter as arguments

    %% Preparation:
    % Arguments:
    if ~exist('n_state','var')
        n_state = 2;
    end
    if ~exist('distribution','var')
        distribution = 'MixtGauss';
    end
    if ~exist('eps_uni','var')
        eps_uni= true;
    end

    % Sizes:
    s_ST = size(set_S{1});
    n_layer = s_ST(2);
    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(set_S{1}{layer}.signal);
    end
    s_image = size(set_S{1}{1}.signal{1});
    n_image = length(set_S);

    % Theta:
    theta = cell(s_ST);

    %% Initialisation values for mean and variance:
    mu = cell(s_ST);
    sigma = cell(s_ST);

    % Go through S to compute mean and variance of the ST:
    % Over all the layers
    for layer=1:n_layer
        % Create the correct number of 'scale':
        mu{layer} = cell(1, n_elmt(1,layer));
        sigma{layer} = cell(1, n_elmt(1,layer));

        % Over all the orientation and scale:
        for scale=1:n_elmt(1,layer)
            mu{layer}{scale} = zeros(s_image);
            % MEAN - Over the 'images' in the set:
            for im=1:n_image
                mu{layer}{scale} = mu{layer}{scale} ...
                    + set_S{im}{layer}.signal{scale};
            end
            mu{layer}{scale} = mu{layer}{scale} / n_image;

            % VARIANCE - Over the 'images' in the set:
            sigma{layer}{scale} = zeros(s_image);
            for im=1:n_image
                sigma{layer}{scale} = sigma{layer}{scale} ...
                    + (set_S{im}{layer}.signal{scale} - mu{layer}{scale}) .^2;
            end
            sigma{layer}{scale} = sigma{layer}{scale} / n_image;
        end
    end

    %% Theta field:
    % Go through S to create the correct data structure:
    % Over all the layers
    for layer=1:n_layer
        % Create the hmm field - a cell of same size as S.signal:
        theta{layer}.proba = cell(1, n_elmt(1,layer));
        theta{layer}.epsilon = cell(1, n_elmt(1,layer));
        theta{layer}.mu = cell(1, n_elmt(1,layer));
        theta{layer}.sigma = cell(1, n_elmt(1,layer));
        theta{layer}.distr = distribution;

        % Over all the orientation and scale:
        for scale=1:n_elmt(1,layer)

            % PROBABILITY:
            % First layer: Multinomial distribution - Tabular CPD :
            if layer == 1
                % P_{S1}(m):
                theta{layer}.proba{scale} = 1/n_state * ones([s_image, n_state]);
            else
                theta{layer}.proba{scale} = 1/n_state * ones([s_image, n_state]);
                %theta{layer}.proba{scale} = zeros([s_im, n_state]);
            end

            % MEANS MU AND VARIANCES SIGMA
            max_mu = max(max(mu{layer}{scale}));
            min_mu = min(min(mu{layer}{scale}));

            for m=1:n_state
                % Mean parameter - as many mean params as states
                theta{layer}.mu{scale}(:,:,m) = mu{layer}{scale} ...
                    + (max_mu-min_mu)/2.*rand(size(mu{layer}{scale})) + min_mu;

                % Variance parameter - as many var params as states
                theta{layer}.sigma{scale}(:,:,m) = sigma{layer}{scale};
            end

            % TRANSITION PROBABILITIES - Tabular CPD
            % For each pixel n_state x n_state posible transitions stored
            % in the 3rd and 4th dimensions of the matrix.
            % The first state (rows) describes the father node, the second
            % state (colums) describe the "studied" child node.
            % e.g.: n_state = 2
            % | (1)L-L  (3)L-H |
            % | (2)H-L  (4)H-H |
            % To access to L-H at layer l and scale j:
            % S{l}.hmm{j}.hidden.epsilon(:,:, 1, 2))
            % Uniform epsilon over father/son:
            if eps_uni
                epsilon = zeros(n_state, n_state);
                
                while any(sum(epsilon,2) ~= ones(n_state,1))
                    % Check that epsilon is summing to 1:
                    % Sometime it seems not to be the case even after
                    % normalisation it is not the case due to rounding error.
                    % (difference of 1e-15)
    %                 while squeeze(sum(theta{layer}.epsilon{scale},4)) ...
    %                             - ones([s_image, n_state]) ~= zeros([s_image, n_state])
                    epsilon = rand(n_state, n_state);
                    rowsum = sum(epsilon,2);
                    epsilon = bsxfun(@rdivide, epsilon, rowsum);
                end
               
                s_eps = [1 1 s_image];
                epsilon = repmat(epsilon, s_eps);
                theta{layer}.epsilon{scale} = permute(epsilon,[3 4 1 2]);
%                 end
            else
                theta{layer}.epsilon{scale} = rand([s_image, n_state, n_state]);
                theta{layer}.epsilon{scale} = ...
                    hmm_prepare_normalEps(theta{layer}.epsilon{scale});
            end
            %+++
            s_check = hmm_Scheck_sum(theta{layer}.epsilon{scale}, ...
                ones(size(sum(theta{layer}.epsilon{scale},4))),...
                'Prep_theta', 'Epsilon', '[1]', layer, scale, verbose);
            if not(s_check)
                a = 0;
            end
            % +++
        end
    end
end


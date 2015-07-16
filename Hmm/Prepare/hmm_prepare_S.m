function [S] = hmm_prepare_S(S, n_state)
% hmm_prepare:  ADD NECESSARY FIELDS TO THE SCATTERING TRANSFORM.
%   Given a scattering transform (cell of structure), this functions 
%   associates to each node the required fields for the HMM modeling.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.  
%
%   - n_state: int
%       Number of states the values of the scattering coefficient will be
%       partitioned into.
%       eg: n_state = 2 --> Small (1) and High (2)
%           n_state = 3 --> Small (1), Medium (2) and High (3)...
%   --------
%   OUTPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.
%       Updated with the appropriate fields:
%       S{layer}.hmm{index}.beta
%       S{layer}.hmm{index}.alpha
%       S{layer}.hmm{index}.state
%       S{layer}.hmm{index}.hmm.parent
%       S{layer}.hmm{index}.hmm.condProb
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - possibility to include priors
%   - use a vargin method
%   - add assert for control
%   - Meta parameter as arguments

    %% Initialization:   
    % Sizes:
    n_layer = length(S);   
    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(S{layer}.signal);
    end
    s_im = size(S{1}.signal{1});
       
    %% Hmm field:
    % Go through S to create the correct data structure:
    % Over all the layers
    for layer=1:n_layer
        % Create the hmm field - a cell of same size as S.signal:
        S{layer}.hmm = cell(size(S{layer}.signal));
        
        % Over all the orientation and scale:
        for scale=1:n_elmt(1,layer)
            
            % CURRENT STATE:
            S{layer}.hmm{scale}.state = randi(n_state, s_im); % Create function sample state
            
            % PARENT:
            S{layer}.hmm{scale}.parent = hmm_find_parent(S, layer, scale);
            
            % CHILDREN:
            S{layer}.hmm{scale}.children = hmm_find_children(S, layer, scale);

            % BROTHERS: 
            S{layer}.hmm{scale}.brothers = hmm_find_brothers(S, layer, scale);

            % CONDITIONAL LIKELIHOODS:
            % beta_i(m)
            S{layer}.hmm{scale}.beta.givenNode = zeros([s_im, n_state]); % delta * ones([s_im, n_state]);
            if layer > 1
            % beta_{i,rho(i)}(m)
                S{layer}.hmm{scale}.beta.givenParents = zeros([s_im, n_state]); % delta * ones([s_im, n_state]);
                % beta_{rho(i)}(m)
                %S{layer}.hmm{scale}.beta.Parents = zeros([s_im, n_state]); %delta * ones([s_im, n_state]);
                % beta_{rho(i)\i}(m)
                S{layer}.hmm{scale}.beta.excludeNode = zeros([s_im, n_state]); %delta * ones([s_im, n_state]); 
            end
            
            % JOINT PROBABILITY FUNCTIONS:
            S{layer}.hmm{scale}.alpha = zeros([s_im, n_state]); %delta * ones([s_im, n_state]);
            
            % CONDITIONAL PROBABILITIES:
            S{layer}.hmm{scale}.condProb.givenVis = zeros([s_im, n_state]); %delta * ones([s_im, n_state]); 
            % Size = [s_im, n_state * n_state] and same indexing for the
            % transition as for the transition probability.
            S{layer}.hmm{scale}.condProb.givenVisAndParents= zeros([s_im, n_state, n_state]);% delta * ones([s_im, n_state * n_state]); % 
        end
    end
end


function [hidStates, dob] = conditional_HIDDEN(S, theta, verbose)           % OK
% conditional_HIDDEN: COMPUTE THE DISTRIBUTION OF THE HIDDEN STATES.
%
%   See "Statistical inference for Hidden MArkov Tree Models and
%   application to wavelet trees" for more details on the algorithm.
%   See algorithm 3 p11.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Structure obtained with the function 'scat' of the 'scatnet' lib.
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib
%       obtained using 'hmm_prepare_theta.
%   - verbose: (optional) bool (default= true)
%       If true then 'hmm_Scheck_sum' displays debugging info.
%
%   --------
%   OUTPUTS:
%   --------
%   - dhidStates: cell(struct)
%       Cell of structures containing the distributions of hidden states.
%   - dob: (optional) struct
%       Debuging OBject, place holder to pass along all variable needed for
%       debugging.
%
%   --------
%   TODO:
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

    n_scale = zeros(1,n_layer);

    % Structure to store the distribution:
    hidStates = cell(1, n_layer);

    for layer=1:n_layer
        n_scale(1,layer) = length(S{layer}.signal);

        % Structure:
        hidStates{layer}.ofHiddenStates = cell(1,n_scale(1,layer));
        % Initialization the matrices:
        for i=1:n_scale(1,layer)
            hidStates{layer}.ofHiddenStates{i} = zeros([s_image n_state]);
        end
    end

    %% Initialisation:
    % Layer 1:
    % P(s_1=k) = Pi_k
    hidStates{1}.ofHiddenStates{1} = theta{1}.proba{1};

    %% Induction:
    % Loop over the layers:
    for layer=2:n_layer
        % Loop over the 'scale' at the layer 'layer':
        for scale=1:n_scale(1,layer)
            % Scale and layer of the father node:
            f_layer = layer-1;
            f_scale = S{layer}.hmm{scale}.parent;

            % [P(s_u=1)...P(s_u=K)] = [P(s_{rho(u)}=1)...P(s_{rho(u)}=K)] *
            % epsilon_{rho(u),u)}
            hidStates{layer}.ofHiddenStates{scale} = ...
                squeeze(sum( ...
                repmat( ...
                hidStates{f_layer}.ofHiddenStates{f_scale}, ...
                1,1,1,n_state) ...
                .* theta{layer}.epsilon{scale} ...
                ,3));

            % NAN, 0 and infinite:
            check_nan = hmm_Scheck_0nan(hidStates{layer}.ofHiddenStates{scale},...
                'cond_HIDDEN', 'distrib.of_hidden_states',...
                layer, scale, verbose);
            % Distrib of Hidden states should sum to [1]:
            check_sum = hmm_Scheck_sum(hidStates{layer}.ofHiddenStates{scale}, ...
                ones(size(sum(hidStates{layer}.ofHiddenStates{scale},3))),...
                'Cond_HIDDEN', 'distrib.of_hidden_states', '[1]', ...
                layer, scale, verbose);
            if not(check_sum)
                a = 0;
            end
            % +++
        end
    end
end


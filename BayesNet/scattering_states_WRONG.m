% WRONG: This function takes the min and the max over a batch of data
% leading to incorrect results

function [minMax_strct, transform ] = scattering_states(transform, n_state)
% scattering_states: TRANSCRIPTS THE SCATTERING COEFFICIENT INTO STATES
%   Given a set of scattering transforms, this functions associate to each
%   coefficient a state among 'n_states' number of possible states based on
%   its value compared to the other coefficient for this element.

%   --------
%   INPUTS:
%   --------
%   - transform: cell
%       Cell of the same length as the number of inputs storing the
%       scattering transform
%   - n_state: (optional) int
%       Number of states the values of the scattering coefficient will be
%       partitioned into.
%       eg: n_state = 2 --> Small (1) and High (2)
%           n_state = 3 --> Small (1), Medium (2) and High (3)
%           ...
%
%   --------
%   OUTPUTS:
%   --------
%   - transform: cell
%       Cell of the same length as the number of inputs storing the
%       scattering transform  updated to contain the 'states'.
%       The states are described by an 'int' where 1 = low and n_state =
%       high

    %% Initialization:
    % Arguments:
    if nargin <2
        n_state = 2;
    end
    
    % Sizes:
    l_batch = length(transform);    
    n_layer = length(transform{1});   
    n_elmt = zeros(1,n_layer);
    
    for l=1:n_layer
        n_elmt(1,l) = length(transform{1}{l}.signal);
    end
    s_st = size(transform{1}{1}.signal{1});
     
    % Max and Min structure:
    minMax_strct = struct();
    minMax_strct.min = {};
    minMax_strct.max = {};
    
    for l=1:n_layer
        minMax_strct.min{l} = {};
        minMax_strct.max{l} = {};
        
        for scl=1:n_elmt(1,l)
            minMax_strct.min{l}{scl} = inf * ones(s_st);
            minMax_strct.max{l}{scl} = -inf * ones(s_st);
        end
    end
    
    % State field in 'transform':
    for im=1:l_batch
        for l=1:n_layer
            transform{im}{l}.state = {};
            for scl=1:n_elmt(1,l)
                transform{im}{l}.state{scl} = zeros(s_st);
            end
        end
    end
        
    %% Minimums and maximums:
    % Over all the images of the batch:
    for im=1:l_batch
        % Over all the scattering layers:
        for l=1:n_layer
            % Over all the orientation and scale:
            for scl=1:n_elmt(1,l)
                minMax_strct.min{l}{scl} = min(minMax_strct.min{l}{scl},...
                                               transform{im}{l}.signal{scl});
                minMax_strct.max{l}{scl} = max(minMax_strct.max{l}{scl},...
                                               transform{im}{l}.signal{scl});       
            end
        end
    end
    
    %% States:
    % Partitioned of the ST coef value space into 'n_states'
    for l=1:n_layer
        % Over all the images of the batch:
        for im=1:l_batch
            % Over all the orientation and scale:
            for scl=1:n_elmt(1,l)
                tmp_interv = (minMax_strct.max{l}{scl} - minMax_strct.min{l}{scl}) ./ n_state;
                
                for s=1:n_state
                    if s == 1
                        transform{im}{l}.state{scl}(transform{im}{l}.signal{scl} <= tmp_interv )= 1;
                    elseif s == n_state
                        transform{im}{l}.state{scl}(transform{im}{l}.signal{scl} >= (n_state-1) .* tmp_interv) = n_state;
                    else
                        transform{im}{l}.state{scl}(...
                            transform{im}{l}.signal{scl} >= (s-1) .* tmp_interv...
                            & ...
                            transform{im}{l}.signal{scl} <= s .*  tmp_interv )...
                                = s;
                    end
                end
            end
        end
    end
end


function [set_S, theta] = hmm_EM(set_S, n_step)
% hmm_EM: COMPUTE THE E"EXPECTATION/MAXIMISATION" ALGORITHM
%   
%   See "Wavelet-based statistical signal processing using hidden markov
%   models" for more details on the algorithm
%
%   --------
%   INPUTS:
%   --------
%   - set_S: cell{cell(struct)}
%       Set of structures obtained with the function 'hmm_EM_E'.      
%   - n_step: int
%       Number of iteration to be done.
%
%   --------
%   OUTPUTS:
%   --------
%   - set_S: cell(struct)
%       Set of structures updated
%   - theta: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib
%       containing the modelization parameters.
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Convergence checking
%   - Argument structure
%   - Input check

    %% Initialization:
    % Arguments:
    distribution = 'MixtGauss';
    n_state = 2;
    
    % Create theta^(0) - Initialisation:
    theta = hmm_prepare_Theta(set_S, n_state, distribution);
    
    % Prepare set_S with the appropriate fields:
    for im=1:length(set_S)
        set_S{im} = hmm_prepare_S(set_S{im}, n_state);
    end
         

    %% EM algo:
    % First step handled sepraratly for timing:
    if n_step >= 1
        tic;
        disp(['--- Step 1/' int2str(n_step) ' ---'])
        set_S = hmm_EM_E(set_S, theta);
        [set_S, theta] = hmm_EM_M(set_S, theta);    
        time = toc;  
        
%         % +++
%         theta{2}.mu{1}(:,:,2)
            
        for l=2:n_step
            % Print:
            disp(['--- Step ' int2str(l) '/' int2str(n_step) ...
                  ' --- expected remaining time: ' ...
                  num2str((n_step-(l-1)) * time) ' s. ---'])

            % E step:
            set_S = hmm_EM_E(set_S, theta);                                 % ----------> Check OK
            % M step:
            [set_S, theta] = hmm_EM_M(set_S, theta);
            
%             % +++
%             theta{2}.mu{1}(:,:,2)
        end
    end
end


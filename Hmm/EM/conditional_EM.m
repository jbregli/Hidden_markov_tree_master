function [ theta, cv_stat, dob] = conditional_EM(set_S, n_step, ...
    n_state, distribution, ...
    eps_uni, verbose, mixing, cv_sens)
% conditional_EM: PERFORM EXPECTATION/MAXIMISATION ON A SET OF IMAGES.
%
%   See "Statistical inference for Hidden MArkov Tree Models and
%   application to wavelet trees" for more details on the algorithm.
%   See algorithm 3 p11.
%
%   --------
%   INPUTS:
%   --------
%   - set_S: cell(cell(struct))
%       Set of structures obtained with the function 'scat' of the
%       'scatnet' lib applied to a set of images.
%   - n_step: int
%       Maximum number of steps for the EM algorithm.
%   - n_state: (optional) int (default: 2)
%       Number of states for the mixture of models.
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
%       Structure with the same organisation as 'S' the 'scatnet' lib
%       containing the modelisation parameters.
%
%   --------
%   TODO:
%   --------
%   - Link the 'model' variable to control which familly of distribution is
%       used for modelisation ('MixtGauss', 'FoldedGauss' ...).

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
    if ~exist('verbose','var')
        verbose = false;
    end
    if ~exist('mixing','var')
        mixing = floor(n_step/10);
    end
    if ~exist('cv_sens','var')
        cv_sens = 1e-3;
    end    
   
    % Variables (2):
    n_image = length(set_S);
    set_hidStates = cell(1,n_image);
    set_proba = cell(1,n_image);

    % Theta:
    [theta, theta_old] = ...
        hmm_prepare_Theta(set_S, n_state, distribution, eps_uni, verbose);
    
    % Sizes (2):
    n_layer = length(theta);
    n_scale = zeros(1,n_layer);
    s_image = size(theta{1}.proba{1}(:,:,1));
    n_state = size(theta{1}.proba{1},3);
    
    % Extreme values test:
    check_strct = cell(1,n_layer);
    
    for layer=1:n_layer
        n_scale(1,layer) = length(theta{layer}.proba);
        check_strct{layer} = cell(1,n_scale(1,layer));
        % Initialization of the matrices:
        for scale=1:n_scale(1,layer)
            check_strct{layer}{scale} = zeros([s_image n_state]);
        end
    end
    
    % Convergence test:
    [cv_ach_strct, cv_ach_bool] = ...
        hmm_conv_test(theta, theta_old, [], 1, mixing, cv_sens, true);        

    % Display:
    fprintf('* EM algorithm: \n');
    reverseStr = '';
    
    % Select a random layer and scale in the scattering transform for 
    % plotting: 
    lay = randi(n_layer);
    scal= randi(n_scale(lay));

    % Select a random pixel in the image:    
    x = randi(s_image(1));
    y = randi(s_image(2));

    %% EM:
    step = 1;

    while (step <= n_step && not(cv_ach_bool))
        % Print remaining steps and times:
        if step == 1
            tic;
            % Display and update:
            msg = sprintf('--- Step %i/%i ---', step, n_step);
            fprintf([reverseStr, msg]);
            if not(verbose)
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        else
            if step ==2
                time = toc;
            end
            % Display and update:
            msg = sprintf('--- Step %i/%i --- Maximum expected remaining time: %.2f s. \r ' ,...
                step, n_step, (n_step-(step-1)) * time);
            fprintf([reverseStr, msg, msg2, msg3, msg4, msg5]);
            if not(verbose)
                reverseStr = repmat(sprintf('\b'), 1,  ...
                    length(msg) + length(msg2) + length(msg3) + length(msg4)...
                    + length(msg5));
            end
        end

        % Old theta for convergence testing:
        if step > 1
            theta_old = theta;
        end
        
        % Expectation:
        % Hidden state probabibility:
        hidStates = conditional_HIDDEN(set_S{1}, theta, verbose);        
        for im=1:n_image
            % UP pass: Compute the betas
            cond_up = conditional_UP(set_S{im}, theta, hidStates, verbose);
            % Down pass: Compute the alphas
            alpha = conditional_DOWN(set_S{im}, theta, hidStates, cond_up, verbose);
            % Conditional probabilities:
            [set_proba{im}, check_strct, s_check] = ...
                conditional_P(set_S{im}, theta, hidStates, ...
                cond_up, alpha, check_strct, cv_ach_strct, verbose);
        end

        % Maximisation:
        % Catch event where the conditional probabilities would be equal to
        % 0.
        if max(s_check) == 0
            theta = ...
                conditional_M(set_S, theta, set_proba, ...
                    eps_uni, cv_ach_strct, theta_old, verbose);
        else
            fprintf('--- Breaking... \n')
            break
        end

        % Convergence testing:
        [cv_ach_strct, cv_ach_bool] = ...
            hmm_conv_test(theta, theta_old, cv_ach_strct, step, mixing,...
                cv_sens); % check_strct);
        
%         % Plot:
%         hmm_plot_distr(set_S, theta, cv_ach_strct, lay, scal, x, y);

        % Display:
        if lay == 1
            msg2 = sprintf('+++ cv_achieved{%i}.proba{%i}=  %s \n', ...
                lay, scal, mat2str(squeeze(cv_ach_strct{lay}.proba{scal}(x,y,:))));
            msg3 = sprintf('+++ theta{%i}.proba{%i} - theta_old{%i}.proba{%i}=  %s \n', ...
                lay, scal, lay, scal, mat2str(squeeze(theta{lay}.proba{scal}(x,y,:)-theta_old{lay}.proba{scal}(x,y,:))));
            msg4 = sprintf('');
            msg5 = sprintf('');
        else
            msg2 = sprintf('+++ cv_achieved{%i}.mu{%i}=  %s \n', ...
                lay, scal, mat2str(squeeze(cv_ach_strct{lay}.mu{scal}(x,y,:))));
            msg3 = sprintf('+++ theta{%i}.mu{%i} - theta_old{%i}.mu{%i}=  %s \n', ...
                lay, scal, lay, scal, mat2str(squeeze(theta{lay}.mu{scal}(x,y,:)-theta_old{lay}.mu{scal}(x,y,:))));
            msg4 = sprintf('+++ cv_achieved{%i}.epsilon{%i}=  %s \n', ...
                lay, scal, mat2str(squeeze(cv_ach_strct{lay}.epsilon{scal}(x,y,:))));
            msg5 = sprintf('+++ theta{%i}.epsilon{%i} - theta_old{%i}.epsilon{%i}=  %s \n', ...
                lay, scal, lay, scal, mat2str(squeeze(theta{lay}.epsilon{scal}(x,y,:)-theta_old{lay}.epsilon{scal}(x,y,:))));
        end
        
        if cv_ach_bool
            break
        end

        % Step iteration
        step=step+1;
    end

    % Convergence:
    if cv_ach_bool
        fprintf('--- Convergence achieved in %i steps. \n', step)
    else
        fprintf('--- Convergence has not yet been achieved after %i steps. \n', ...
            step-1)
    end

    % Statistics on the number of converged pixels:
    cv_stat = hmm_conv_stat(cv_ach_strct);
    
    % +++ Debuging object:
    dob.set_distrib = set_hidStates;
    dob.cond_up = cond_up;
    dob.alpha = alpha;
    dob.set_proba = set_proba;
end


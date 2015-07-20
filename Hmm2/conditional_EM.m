function [ theta, dob] = conditional_EM(set_S, n_step, ...
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
        disp('nop')
        eps_uni= true;
    end
    if ~exist('verbose','var')
        verbose = true;
    end
    if ~exist('mixing','var')
        mixing = floor(n_step/10);
    end
    if ~exist('cv_sens','var')
        cv_sens = 1e-4;
    end

    % Variables:
    n_image = length(set_S);
    set_hidStates = cell(1,n_image);
    set_proba = cell(1,n_image);

    % Theta:
    [theta, theta_old] = ...
        hmm_prepare_Theta(set_S, n_state, distribution, eps_uni, verbose);
    
    % Convergence test:
    [cv_ach_strct, cv_ach_bool] = ...
        hmm_conv_test(theta, theta_old, 1, 'init', 'true');

    % Display:
    fprintf('* EM algorithm: \n');
    reverseStr = '';

    %% EM:
    step = 1;

    while (step <= n_step && not(cv_ach_bool))
        % Print remaining steps and times:
        if step == 1
            tic;
            msg = sprintf('--- Step %i/%i ---', step, n_step);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        else
            if step ==2
                time = toc;
            end
            msg = sprintf('--- Step %i/%i --- Maximum expected remaining time: %.2f s. \r ' ,...
                step, n_step, (n_step-(step-1)) * time);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));

        end

        % Old theta for convergence testing:
        if step > 1
            theta_old = theta;
        end

        % s_check matrix:
        s_check = ones(1,n_image);

        % Expectation:
        for im=1:n_image
            set_hidStates{im} = conditional_HIDDEN(set_S{im}, theta, verbose);
            cond_up = conditional_UP(set_S{im}, theta, set_hidStates{im}, verbose);
            alpha = conditional_DOWN(set_S{im}, theta, set_hidStates{im}, cond_up, verbose);
            [set_proba{im}, s_check(1,im)] = ...
                conditional_P(set_S{im}, theta, set_hidStates{im}, cond_up, alpha, verbose);
        end

        % Maximisation:
        if min(s_check) == 1
            theta = conditional_M(set_S, theta, set_hidStates, set_proba, ...
                eps_uni, cv_ach_strct, theta_old, verbose);
        else
            fprintf('--- Breaking... \n')
            break
        end

        % Convergence testing:
        [cv_ach_strct, cv_ach_bool] = hmm_conv_test(theta, theta_old, step, mixing, ...
            cv_sens);

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

    % +++ Debuging object:
    dob.set_distrib = set_hidStates;
    dob.cond_up = cond_up;
    dob.alpha = alpha;
    dob.set_proba = set_proba;
end


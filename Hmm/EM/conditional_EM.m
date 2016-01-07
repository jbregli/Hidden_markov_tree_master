function [ theta, cv_stat, dob] = conditional_EM(set_S, ...
    EM_metaparameters, options)

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
%   - EM_metaparameters: (optional) struct
%        Each fields is a meta parameter for the algo. 
%           - .n_step: (optional) int (default: 100)
%               Maximum number of steps for the EM algorithm.
%           - .n_state: (optional) int (default: 2)
%               Number of states for the mixture of models.
%           - .distribution: (optional) str (default: MixtGauss)
%               Distribution to be used for modeling.
%               Options: MixtGauss
%           - .eps_uni: (optional) bool (default: false)
%               If set to true epsilon is considered to be uniform accross 
%               the image (ie: same probability transition from father to
%               child for all the pixel of a same image.
%           - .mixing: (optional) int (default: 10)
%               Number of step beforetesting convergence.
%           - .cv_sens: (optional) float (default: 1e-5)
%               Minimum increase between two steps to consider convergence.
%           - .cv_steps: (optional) int (default: 7)
%               Number of consecutive steps where the increase has to be 
%               lower than 'cv_sens' to consider convergence.
%           - .cv_ratio: (optional) float[0,1] (default: 0.98) 
%                Ratio of converged pixel before breaking.
%   - options: (optional) struct
%        Each fields is an option for the algo. 
%           - .verbose: (optional) bool (default= false)
%                   If true then 'hmm_Scheck_sum' displays debugging info.
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
    if nargin < 2
        EM_metaparameters = struct();
    end
    if nargin < 3
        options = struct();
    end        
    
    % Fields of the input structure:
    EM_metaparameters = isfield_EMmeta(EM_metaparameters);
    options = isfield_options(options);

    
    % Variables (2):
    n_image = length(set_S);
    set_hidStates = cell(1,n_image);
    set_proba = cell(1,n_image);

    % Theta:
    [theta, theta_old] = ...
        hmm_prepare_Theta(set_S, EM_metaparameters.n_state, EM_metaparameters.distribution, EM_metaparameters.eps_uni, options.verbose);
    
    % Sizes (2):
    n_layer = length(theta);
    n_scale = zeros(1,n_layer);
    s_image = size(theta{1}.proba{1}(:,:,1));
    EM_metaparameters.n_state = size(theta{1}.proba{1},3);
    
    % Extreme values test:
    check_strct = cell(1,n_layer);
    
    for layer=1:n_layer
        % Number of scale: per layer:
        n_scale(1,layer) = length(theta{layer}.proba);
        
        check_strct{layer}.ofNode = cell(1,n_scale(1,layer));
        check_strct{layer}.ofNodeAndParent = cell(1,n_scale(1,layer));
        
        % Initialization of the matrices:
        for scale=1:n_scale(1,layer)
            check_strct{layer}.ofNode{scale} = zeros([s_image EM_metaparameters.n_state]);
            check_strct{layer}.ofNodeAndParent{scale} = ...
                zeros([s_image EM_metaparameters.n_state EM_metaparameters.n_state]);
        end
    end
    
    % Convergence test:
    [cv_status, cv_ach_bool] = ...
        hmm_conv_test(theta, theta_old, [], 1, ...
        EM_metaparameters.mixing, EM_metaparameters.cv_sens, ...
        EM_metaparameters.cv_steps, true);        

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

    while (step <= EM_metaparameters.n_step && not(cv_ach_bool))
        % Print remaining steps and times:
        if step == 1
            tic;
            % Display and update:
            msg = sprintf('--- Step %i/%i ---', step, EM_metaparameters.n_step);
            fprintf([reverseStr, msg]);
            if not(options.verbose)
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        else
            if step == 2
                time = toc;
            end
            % Display and update:
            if EM_metaparameters.n_step == inf
                msg = sprintf('--- Step %i --- Single step time: %.2f s. \r ' ,...
                    step, time);
            else
                msg = sprintf('--- Step %i/%i --- Maximum expected remaining time: %.2f s. \r ' ,...
                    step, EM_metaparameters.n_step, (EM_metaparameters.n_step-(step-1)) * time);
            end
            
            fprintf([reverseStr, msg, msg2, msg3, msg4, msg5, msg6]);
            if not(options.verbose)
                reverseStr = repmat(sprintf('\b'), 1,  ...
                    length(msg) + length(msg2) + length(msg3) + length(msg4)...
                    + length(msg5) + length(msg6));
            end
        end

        % Old theta for convergence testing:
        if step > 1
            theta_old = theta;
        end
        
        % Expectation:
        % Hidden state probabibility:
        hidStates = conditional_HIDDEN(set_S{1}, theta, options.verbose);        
        for im=randperm(n_image)
            % UP pass: Compute the betas
            cond_up = conditional_UP(set_S{im}, theta, hidStates, options.verbose);
            % Down pass: Compute the alphas
            alpha = conditional_DOWN(set_S{im}, theta, hidStates, cond_up, options.verbose);
            % Conditional probabilities:
            [set_proba{im}, check_strct, ~] = ...
                conditional_P(set_S{im}, theta, hidStates, ...
                cond_up, alpha, check_strct, options.verbose);
        end

        % Maximisation:
        % Catch event where the conditional probabilities would be equal to
        theta = ...
            conditional_M(set_S, theta, set_proba, ...
                EM_metaparameters.eps_uni, cv_status, theta_old, ...
                check_strct, options.verbose);

        % Check theta values:
        check_bool = hmm_check_theta(theta);
        if check_bool
            fprintf('--- Breaking - theta test... \n')
            theta = theta_old;
            break
        end
                
        % Convergence testing:
        [cv_status, cv_ach_bool] = ...
            hmm_conv_test(theta, theta_old, cv_status, step, ...
                EM_metaparameters.mixing, EM_metaparameters.cv_sens, ...
                EM_metaparameters.cv_steps, false);
        
        % +++ 
        % Statistics on the number of converged pixels:
        cv_stat = hmm_conv_stat(cv_status);
        
        % Display:
        if lay == 1
            msg2 = sprintf('+++ cv_status{%i}.proba{%i}=  %s \n', ...
                lay, scal, mat2str(squeeze(cv_status.params{lay}.proba{scal}(x,y,:))));
            msg3 = sprintf('+++ theta{%i}.proba{%i} - theta_old{%i}.proba{%i}=  %s \n', ...
                lay, scal, lay, scal, mat2str(squeeze(theta{lay}.proba{scal}(x,y,:)-theta_old{lay}.proba{scal}(x,y,:))));
            msg4 = sprintf('');
            msg5 = sprintf('');
            msg6 = sprintf('');
        else
            msg2 = sprintf('+++ cv_status{%i}.mu{%i}=  %s \n', ...
                lay, scal, mat2str(squeeze(cv_status.params{lay}.mu{scal}(x,y,:))));
            msg3 = sprintf('+++ theta{%i}.mu{%i} - theta_old{%i}.mu{%i}=  %s \n', ...
                lay, scal, lay, scal, mat2str(squeeze(theta{lay}.mu{scal}(x,y,:)-theta_old{lay}.mu{scal}(x,y,:))));
            msg4 = sprintf('+++ cv_status{%i}.epsilon{%i}=  %s \n', ...
                lay, scal, mat2str(squeeze(cv_status.params{lay}.epsilon{scal}(x,y,:))));
            msg5 = sprintf('+++ theta{%i}.epsilon{%i} - theta_old{%i}.epsilon{%i}=  %s \n', ...
                lay, scal, lay, scal, mat2str(squeeze(theta{lay}.epsilon{scal}(x,y,:)-theta_old{lay}.epsilon{scal}(x,y,:))));
        end
        
        msg6 = sprintf('+++ cv_stat.overallCv = %.5f \n', cv_stat.overallCv);
        
        if cv_stat.overallCv > EM_metaparameters.cv_ratio
            cv_ach_bool = true;
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
        if EM_metaparameters.rerun == true
            
            fprintf('RE-RUN %i \n', EM_metaparameters.rerun_count)
            EM_metaparameters.rerun_count = ...
                EM_metaparameters.rerun_count +1 ;
            [theta, ~, ~] = conditional_EM(set_S, ...
                EM_metaparameters, options);
        end
    end

    % Statistics on the number of converged pixels:
    cv_stat = hmm_conv_stat(cv_status);
    
    % +++ Debuging object:
    dob.set_distrib = set_hidStates;
    dob.cond_up = cond_up;
    dob.alpha = alpha;
    dob.set_proba = set_proba;
    dob.cv_status = cv_status;
end


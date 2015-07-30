%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realizes a convergence test of the EM algorithm.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We create a simulated probabilistic tree:
% Initial state distribution:

%clear all
%close all

%% Initialization:
real_data = true;

% Size of the simulated images:
s_im = [10 10];
n_state = 2;

% Number of "images" in the set;
n_image = 50;

% Number of optimization step:
n_step = inf;

% Model distribution:
distribution = 'MixtGauss';
% Epsilon uniform over the pixels of a father/son transition
eps_uni= false;
% Display error messages:
verbose = true;
% Sensibility f the convergence test:
cv_sens = 1e-6;

%% Creating the proper tree structure:
% The scattering coefficients are not used yet. 
% Later in the code images are sampled from ground truth.
% REAL DATA
if real_data
%     directory = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';
%     label = 'corduroy/'; 
%     path_to_set = fullfile(directory, label);
    path_to_set = {'square', n_image, s_im};
else
    path_to_set = {'square', n_image, s_im};
end

% Parameters:
filt_opt.J = 2; % scales
filt_opt.L = 2; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;
% filt_opt = struct();
% scat_opt = struct();

% ST:
set_S = scat_class(path_to_set, filt_opt, scat_opt);
S = set_S{1};

% REAL DATA:
if real_data
    for im=1:length(set_S)
        set_S{im} = hmm_prepare_S(set_S{im}, n_state);
    end
else
    %% Theta for simultation:
    % Each node is a mixture of 2 gaussians:
    mu_L = 1; sigma_L = .3;
    mu_H = 5; sigma_H = .5;
    
    if n_state == 2
        dist_param = {{mu_L sigma_L} {mu_H sigma_H}};
        
        % Transition proba:
        % | (1)L-L  (3)L-H |
        % | (2)H-L  (4)H-H |
        epsilon = [0.8 0.1; 0.2 0.9]; 
        
        % Root node proba:
        rn_prob = [0.5 0.5];
        
    elseif n_state == 3
        mu_M = 2; sigma_M = .2;
        dist_param = {{mu_L sigma_L} {mu_M sigma_M} {mu_H sigma_H}};

        % Transition proba:
        epsilon = [0.7 0.2 0.1; 0.1 0.7 0.2; 0.2 0.1 0.7];

        % Root node proba:
        rn_prob = [0.3 0.3 0.4];
    else
        fprintf('This script can only handle simultations with 2 or 3 gaussian distribution mixed')
    end

    %% Reformat the tree:
    for l=1:length(S)
        for s=1:length(S{l}.signal)
            S{l}.signal{s} = zeros(s_im);
        end
    end
    S = hmm_prepare_S(S, n_state);

    %% Sample the images:
    theta_gt = GT_theta(S, dist_param, epsilon, rn_prob);
    set_S = GT_set(S, theta_gt, n_image);
end

%% Hmm model:
[theta_est, cv_stat, dob] = conditional_EM(set_S, n_step, n_state, distribution, ...
                              eps_uni, verbose, 10, cv_sens);

if real_data == false
    error_LW = hmm_error(theta_est, theta_gt);
end
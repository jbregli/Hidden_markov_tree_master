%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realizes a convergence test of the EM algorithm.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We create a simulated probabilistic tree:
% Initial state distribution:

%clear all
%close all

%% Initialization:
real_data = false;

% Size of the simulated images:
s_im = [10 10];
n_state = 2;

% Number of "images" in the set;
n_image = 20;

%Number of optimization step:
n_step = 100;

% Model distribution:
distribution = 'MixtGauss';
% Epsilon uniform over the pixels of a father/son transition
eps_uni= false;
% Display error messages:
verbose = true;

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
        set_S{im} = hmm_prepare_S( set_S{im}, n_state);
    end
else
    %% Theta for simultation:
    % Each node is a mixture of 2 gaussians:
    mu_L = 1; sigma_L = .1;
    mu_H = 5; sigma_H = .1;
    dist_param = {{mu_L sigma_L} {mu_H sigma_H}};

    % Transition proba:
    % | (1)L-L  (3)L-H |
    % | (2)H-L  (4)H-H |
    epsilon = [0.7 0.1; 0.3 0.9];

    % Root node proba:
    rn_prob = [0.3 0.7];

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

% REAL DATA:
% for im=1:length(set_S)
%     set_S{im} = hmm_prepare_S( set_S{im}, n_state);
% end

%% Hmm model:
[theta, dob] = conditional_EM(set_S, n_step, n_state, distribution, ...
                              eps_uni, verbose);
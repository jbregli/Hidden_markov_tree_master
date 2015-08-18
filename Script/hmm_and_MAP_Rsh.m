%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realizes a convergence test of the EM algorithm.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We create a simulated probabilistic tree:
% Initial state distribution:

clear all
close all

%% Initialization:
% Size of the simulated images:
s_im = [10 10];
n_state = 2;

% Number of "images" in the set;
n_image = 100;

% Number of optimization step:
n_step = 100;

% Model distribution:
distribution = 'MixtGauss';
% Epsilon uniform over the pixels of a father/son transition
eps_uni= false;
% Display error messages:
verbose = false;
% Sensibility f the convergence test:
cv_sens = 1e-6;

%% Creating the proper tree structure:
% The scattering coefficients are not used yet. 
% Later in the code images are sampled from ground truth.
% REAL DATA
path_to_set = {'square', n_image, s_im};

% Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;
% filt_opt = struct();
% scat_opt = struct();

% ST:
set_S = scat_class(path_to_set, filt_opt, scat_opt);
S = set_S{1};

% Prepare the scattering structure for HMM:
for im=1:length(set_S)
    set_S{im} = hmm_prepare_S(set_S{im}, n_state);
end

%% Hmm model:
[theta_est, cv_stat, dob] = conditional_EM(set_S, n_step, n_state, distribution, ...
                              eps_uni, verbose, 10, cv_sens);

%% Error of modelisation:
% Need to find something for real data
% if real_data == false
%     error_LW = hmm_error(theta_est, theta_gt);
% end

%% MAP:
% SQUARE:
% Creating a new square image:
path_to_square = {'square', 1, s_im};
% ST:
square_S = scat_class(path_to_square, filt_opt, scat_opt);
S_square = square_S{1};

% Prepare the scattering structure for HMM:
for im=1:length(square_S)
    square_S{im} = hmm_prepare_S(square_S{im}, n_state);
end

[P_hat_square, H_tree_square] = hmm_MAP(square_S{1}, theta_est, verbose);

% CIRCLE:
path_to_circle = {'circle', 1, s_im};
% ST:
circle_S = scat_class(path_to_circle, filt_opt, scat_opt);
S_circle = circle_S{1};

% Prepare the scattering structure for HMM:
for im=1:length(circle_S)
    circle_S{im} = hmm_prepare_S(circle_S{im}, n_state);
end

[P_hat_circle, H_tree_circle] = hmm_MAP(circle_S{1}, theta_est, verbose);
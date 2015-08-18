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
directory = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';
label = 'corduroy/'; 
path_to_corduroy = fullfile(directory, label);

% Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;
% filt_opt = struct();
% scat_opt = struct();

% ST:
corduroy_S = scat_class(path_to_corduroy, filt_opt, scat_opt);
S = corduroy_S{1};

% Prepare the scattering structure for HMM:
for im=1:length(corduroy_S)
    corduroy_S{im} = hmm_prepare_S(corduroy_S{im}, n_state);
end

% Remove the last element to be kept for testing:
corduroy_T = corduroy_S{end};
corduroy_S = corduroy_S(1:end-1);
 
%% Hmm model:
[theta_est, cv_stat, dob] = conditional_EM(corduroy_S, n_step, n_state, distribution, ...
                              eps_uni, verbose, 10, cv_sens);

%% Error of modelisation:
% Need to find something for real data
% if real_data == false
%     error_LW = hmm_error(theta_est, theta_gt);
% end

%% MAP:
% Corduroy
[P_hat_cord, H_tree_cord] = hmm_MAP(corduroy_T, theta_est, verbose);

% CIRCLE:
label = 'orange_peel'; 
path_to_orange = fullfile(directory, label);

% ST:
orange_T = scat_class(path_to_orange, filt_opt, scat_opt,1);

% Prepare the scattering structure for HMM:
for im=1:length(orange_T)
    orange_T{im} = hmm_prepare_S(orange_T{im}, n_state);
end

[P_hat_oran, H_tree_oran] = hmm_MAP(orange_T{1}, theta_est, verbose);

msg_cord = sprintf('P(corduroy = corduroy) = %.5f \r', ...
    mean(mean(P_hat_cord)));
msg_oran = sprintf('P(orange = corduroy) = %.5f \r', ...
    mean(mean(P_hat_oran)));

fprintf(msg_cord, msg_oran)

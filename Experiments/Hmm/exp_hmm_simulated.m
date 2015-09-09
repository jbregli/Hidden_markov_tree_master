%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test of the STHMT on genereted    %
% trees with known parameters (intial state distrib, transition proba and %
% gaussian means and variances.                                           %
% The SCHMT is trained on one set of parameters.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:
% Size of the simulated images:
s_im = [10 10];
n_state = 2;

% Number of "images" in the set;
n_image = 100;

% Number of optimization step:
n_step = 200;

% Model distribution:
distribution = 'MixtGauss';
% Epsilon uniform over the pixels of a father/son transition
eps_uni= true;
% Display error messages:
verbose = false;
% Sensibility f the convergence test:
cv_sens = 1e-6;

%% Creating the proper tree structure:
% The scattering coefficients are not used yet. 
% Later in the code images are sampled from ground truth.
path_to_set = {'square', n_image, s_im};

% Parameters:
filt_opt.J = 4; % scales
filt_opt.L = 4; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;
% filt_opt = struct();
% scat_opt = struct();

% ST:
set_S = ST_class(path_to_set, filt_opt, scat_opt);
S = set_S{1};

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
theta_gt = create_theta_groundT(S, dist_param, epsilon, rn_prob);
set_S = create_set_simul(S, theta_gt, n_image);


%% Hmm model:
[theta_est, cv_stat, dob] = conditional_EM(set_S, n_step, n_state, distribution, ...
                              eps_uni, verbose, 10, cv_sens);


%% Modelisation error:                          
error_LW = hmm_error(theta_est, theta_gt);

%% MAP:
% Create a new image from the same distribution:
test_S_T = create_set_simul(S, theta_gt, 1);

[P_hat_T, H_tree_T] = hmm_MAP(test_S_T{1}, theta_est, verbose);

% Create a new image from another distribution:
mu_L2 = 2; sigma_L2 = .1;
mu_H2 = 10; sigma_H2 = .5;
if n_state == 2
    dist_param2 = {{mu_L2 sigma_L2} {mu_H2 sigma_H2}};

    % Transition proba:
    epsilon2 = [0.5 0.5; 0.5 0.5]; 

    % Root node proba:
    rn_prob2 = [0.5 0.5];

elseif n_state == 3
    mu_M2 = 5; sigma_M2 = .2;
    dist_param2 = {{mu_L2 sigma_L2} {mu_M2 sigma_M2} {mu_H sigma_H2}};

    % Transition proba:
    epsilon2 = [0.4 0.3 0.3; 0.3 0.5 0.3; 0.3 0.3 0.4];

    % Root node proba:
    rn_prob2 = [0.3 0.3 0.4];
else
    fprintf('This script can only handle simultations with 2 or 3 gaussian distribution mixed')
end

theta_gt2 = create_theta_groundT(S, dist_param2, epsilon2, rn_prob2);
test_S_F = create_set_simul(S, theta_gt2, 1);

[P_hat_F, H_tree_F] = hmm_MAP(test_S_F{1}, theta_est, verbose);

msg_T = sprintf('P(M1 = M1) = %.5f \r', ...
    mean(mean(P_hat_T)));
msg_F = sprintf('P(M2 = M1) = %.5f \r', ...
    mean(mean(P_hat_F)));

fprintf(msg_T);
fprintf(msg_F);
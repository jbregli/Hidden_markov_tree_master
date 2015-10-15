%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test of the STHMT on genereted    %
% shapes (square and circle).                                             %
% The SCHMT is trained to model squares.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Initialization:
% Size of the simulated images:
s_im = [30 30];

% Number of states:
n_state = 2;

% Number of "images" in the training set;
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

% Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 6; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;


%% CREATE THE TREE STRUCTURE
% The scattering coefficients are not used yet. 
% Later in the code images are sampled from ground truth.
path_to_set = {'square', 1, s_im};

% ST:
tmp_set = ST_class(path_to_set, filt_opt, scat_opt);
S = tmp_set{1};

% Reformat the tree:
for l=1:length(S)
    for s=1:length(S{l}.signal)
        S{l}.signal{s} = zeros(s_im);
    end
end
S = hmm_prepare_S(S, n_state);


%% CLASS 1 - DISTRIBUTION 1 -PARAMETERS & SAMPLING:
% Theta for simultation:
% Each node is a mixture of 2 gaussians:
mu_L = 1; sigma_L = .3;
mu_H = 5; sigma_H = .5;

dist_param1 = {{mu_L sigma_L} {mu_H sigma_H}};

% Transition proba:
% | (1)L-L  (3)L-H |
% | (2)H-L  (4)H-H |
epsilon = [0.8 0.1; 0.2 0.9]; 

% Root node proba:
rn_prob = [0.5 0.5];

% Sample the images:
theta_gt_D1 = create_theta_groundT(S, dist_param1, epsilon, rn_prob);
set_S_D1 = create_set_simul(S, theta_gt_D1, n_image);


%% CLASS 2 - DISTRIBUTION 2 - PARAMETERS& SAMPLING:
% Theta for simultation:
% Create a new image from another distribution:
mu_L2 = 2; sigma_L2 = .1;
mu_H2 = 10; sigma_H2 = .5;

dist_param2 = {{mu_L2 sigma_L2} {mu_H2 sigma_H2}};

% Transition proba:
epsilon2 = [0.5 0.5; 0.5 0.5]; 

% Root node proba:
rn_prob2 = [0.5 0.5];

% Sample the images:
theta_gt_D2 = create_theta_groundT(S, dist_param2, epsilon, rn_prob);
set_S_D2 = create_set_simul(S, theta_gt_D1, n_image);

%% CLASS 1 - DISTIRBUTION 1 - TRAINING: 
fprintf('------ TRAINING DISTRIBUTION 1 ------ \n')

[theta_est_D1, cv_stat_D1, ~] = ...
    conditional_EM(set_S_D1, n_step, n_state, distribution, eps_uni, ...
        verbose, 10, cv_sens);

%% CLASS 2 - DISTIBUTION 2 - TRAINING: 
fprintf('------ TRAINING DISTRIBUTION 2 ------ \n')

[theta_est_D2, cv_stat_D2,  ~] = ...
    conditional_EM(set_S_D2, n_step, n_state, distribution, eps_uni, ...
        verbose, 10, cv_sens);    
    
    
%% MAP - CLASSIFICATION SCORE:
fprintf('------ TESTING ------ \n')
n_test = 50;
score_D1= 0;
score_D2 = 0;


P_hat_D1_D1 = cell(1,n_test);
P_hat_D1_D2 = cell(1,n_test);
P_hat_D2_D1 = cell(1,n_test);
P_hat_D2_D2 = cell(1,n_test);

% Creating the test set:
D1_S = create_set_simul(S, theta_gt_D1, n_test);
D2_S = create_set_simul(S, theta_gt_D2, n_test);

for im=1:n_test
    % MAP model = D1:
    [tmp_P_hat_D1_D1, ~] = ...
        hmm_MAP(D1_S{im}, theta_est_D1, false);
    [tmp_P_hat_D1_D2, ~] = ...
        hmm_MAP(D1_S{im}, theta_est_D2, false);
    
    P_hat_D1_D1{im} = mean(mean(tmp_P_hat_D1_D1));
    P_hat_D1_D2{im} = mean(mean(tmp_P_hat_D1_D2));

    if max(P_hat_D1_D2{im}, P_hat_D1_D1{im}) == P_hat_D1_D1{im}
        score_D1 = score_D1 + 1/n_test;
    end
    
    % MAP model = circle:
    [tmp_P_hat_D2_D1, ~] = ...
        hmm_MAP(D2_S{im}, theta_est_D1, false);
    [tmp_P_hat_D2_D2, ~] = ...
        hmm_MAP(D2_S{im}, theta_est_D2, false);
    
    P_hat_D2_D1{im} = mean(mean(tmp_P_hat_D2_D1));
    P_hat_D2_D2{im} = mean(mean(tmp_P_hat_D2_D2));

    if max(P_hat_D2_D2{im}, P_hat_D2_D1{im}) == P_hat_D2_D2{im}
        score_D2 = score_D2 + 1/n_test;
    end
end

fprintf('The total classification score is %.4f. \n', ...
    (score_D2 + score_D1)/2)
fprintf('The classification score for distribution 1 is %.4f. \n', ...
    score_D1)
fprintf('The total classification score for distribution 2  is %.4f. \n', ...
    score_D2)
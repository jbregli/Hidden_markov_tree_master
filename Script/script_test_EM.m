%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realizes a convergence test of the EM algorithm.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We create a simulated probabilistic tree:
% Initial state distribution:

%clear all
%close all

%% Initialization:
% if not(exist('add_path','var'))
%     add_path = 1;
%     addpath_hmt
% end

% Size of the simulated images:
s_im = [10 10];
n_state = 2;

% Number of "images" in the set;
n_image = 300;

%Number of optimization step:
n_step = 40;

%% Creating the proper tree structure:
% The scattering coefficients are not used yet. 
% Later in the code images are sampled from ground truth.
path_to_set = {'square', 1, s_im}; %fullfile(directory, label);

% Parameters:
filt_opt.J = 2; % scales
filt_opt.L = 1; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 1;
% filt_opt = struct();
% scat_opt = struct();

% ST:
set_S = scat_class(path_to_set, filt_opt, scat_opt);
S = set_S{1};

%% Theta for simultation:
% Each node is a mixture of 2 gaussians:
mu_L = 0; sigma_L = 0.1;
mu_H = 0; sigma_H = 1;
dist_param = {{mu_L sigma_L} {mu_H sigma_H}};

% Transition proba:
% | (1)L-L  (3)L-H |
% | (2)H-L  (4)H-H |
epsilon = [0.6 0.1; 0.4 0.9];

% Root node proba:
rn_prob = [0.4 0.6];

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

%% Hmm model:
[set_S, theta] = hmm_EM(set_S, n_step);



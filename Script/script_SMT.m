%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realizes the scattering transform of several images of a    %
% same class. Then it fit a markov tree to the representation created.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possible improvements:
%   - Check if the parameters of the ST have been changed when skipping
%     recomputing

%clear all
close all

%% Initialization:
if not(exist('add_path','var'))
    add_path = 1;
    addpath_hmt
end

% Dataset:
%directory = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';

%% Images from a class:
label = 'corduroy/'; 
% aluminium_foil corduroy cracker orange_peel sandpaper styrofoam 
% brown_bread cotton linen sponge

path_to_set = {'square', 20, [100, 100]}; %fullfile(directory, label);



%% Scattering transform:
% Parameters:
filt_opt.J = 4; % scales
filt_opt.L = 2; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 3;
% filt_opt = struct();
% scat_opt = struct();

% ST - computed only if tansform doesn't exist
% exst = exist('transform_BN','var')
if not(exist('set_S','var'))
    set_S = scat_class(path_to_set, filt_opt, scat_opt);
end

%% Create the HMT:
n_step = 50;

[set_S, theta] = hmm_EM(set_S, n_step);

% % for the first image:
% distribution = 'MixtGauss';
% n_state = 2;
% 
% S = hmm_prepare(set_S{1}, n_state, distribution);
% 
% % Test up - TB applied to all the image zithin the set
% S = hmm_EM_E_up(S);
% % Test down:
% S = hmm_EM_E_down(S);
% % % Test conditional probabilities:
% S= hmm_EM_E_CP(S);
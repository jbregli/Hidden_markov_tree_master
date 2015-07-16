%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script realizes the scattering transform of several images of a    %
% same class. Then it plots the histogram of the ST coefficient           %
% distributions.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:
addpath_hmt

% Dataset:
directory = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';

% Plot info:
layer = 2;
index = 10;
pixel = 250;

%% Images from a class:
label = 'corduroy/'; 
% aluminium_foil corduroy cracker orange_peel sandpaper styrofoam 
% brown_bread cotton linen sponge

path_to_set = fullfile(directory, label);

%% Scattering transform:
% Parameters:
filt_opt.J = 4; % scales
filt_opt.L = 4; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 3;
% filt_opt = struct();
% scat_opt = struct();

% ST:
transform = scat_class(path_to_set, filt_opt, scat_opt);

%% Display:
data = plot_ST_histo(transform, 2, 10, 500);
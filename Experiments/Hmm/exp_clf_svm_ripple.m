%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test of the STHMT on genereted    %
% cropped sonar images from MUSSEL AREA C dataset.                        %
% The SCHMT is trained to model squares.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Initialization:
% Data path:
dir_training = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Training/';

%% CLASSES - RIPPLE 0 - SEABED 1: 
label_ripple = 'Ripple/'; 
path_to_training_ripple = fullfile(dir_training, label_ripple);
label_seabed = 'Mix'; 
path_to_training_seabed = fullfile(dir_training, label_seabed);

% List all the images:
[data_ripple, group_ripple ] = load_svm(path_to_training_ripple, 0);
[data_seabed, group_seabed ] = load_svm(path_to_training_seabed, 1);

data_training = vertcat(data_ripple, data_seabed);
group_training = vertcat(group_ripple, group_seabed);

%% TRAINING:
fprintf('------ TRAINING ------ \n')

SVMStruct = fitcsvm(data_training, group_training);

%% TESTING - CLASSIFICATION SCORE:
fprintf('------ TESTING ------ \n')

dir_test = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Test/';

path_to_test_ripple = fullfile(dir_test, label_ripple);
path_to_test_seabed = fullfile(dir_test, label_seabed);

[data_ripple, group_ripple ] = load_svm(path_to_test_ripple, 0);
[data_seabed, group_seabed ] = load_svm(path_to_test_seabed, 1);

data_test = vertcat(data_ripple, data_seabed);
group_test = vertcat(group_ripple, group_seabed);

label = predict(SVMStruct, data_test);

score = sum(abs(label-group_test)) / length(group_test);
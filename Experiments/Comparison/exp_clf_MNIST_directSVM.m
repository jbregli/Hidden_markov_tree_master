%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test on few classes from the      %
% Kylberg Non rotated texture dataset.                                    %
%                                                                         %
% BEST PARAMETER SO FAR: CS=0.8                                           %
% n_image = 0; n_state = 2; n_step = 100; eps_uni= false; cv_sens = 1e-5; %
% filt_opt.J = 5; filt_opt.L = 3; filt_opt.filter_type = 'morlet';        %
% scat_opt.oversampling = 2; scat_opt.M = 2;                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
close all


%% ===== Initialization: =====
n_training = 5; % Dataset training + testing has 160 images per class 
n_testing = 20;

% SVM parameters:
svm.model = '-g 0.0018 -c 3';


rdm_spl = randsample(1:5000, 5000);

rdm_training = rdm_spl(1:n_training);
rdm_testing = rdm_spl(n_training+1:n_training+ n_testing); % rdm_spl(n_training+1:end);

%% ===== PREPARE DATA: =====
% Data path:
im_format = 'ubyte';
dir_data = '/home/jeanbaptiste/Datasets/Mnist/';
addpath(dir_data);

% Loard MNIST:
tmp_mnist = loadMNISTImages('train-images-idx3-ubyte');
tmp_label = loadMNISTLabels('train-labels-idx1-ubyte');


% Select n_train training example per label
features = []; test_features = [];
labels = []; test_labels = [];
for label=0:n_label-1
    tmp_feature = tmp_mnist(:,tmp_label==label);
    tmp_feature_train = tmp_feature(:, rdm_training);
    tmp_feature_test = tmp_feature(:, rdm_testing);
    
    features = vertcat(features, tmp_feature_train);
    labels = vertcat(labels, repmat(label,n_training,1));
    
    test_features = horzcat(test_features, tmp_feature_test);
    test_labels = vertcat(test_labels, repmat(label,n_testing,1));
end
features  = features';
test_features = test_features';

%% ===== SVM TRAINING: =====
fprintf('------ TRAINING ------ \n')

% Shuffle elements
shuffler = randperm(size(features,1));
features = features(shuffler,:);
labels = labels(shuffler,:);

% Training:
options = {};
svm_params = svm_cross_validation(labels, double(features), options);
svm_model = svmtrain(labels, double(features), svm_params);


%% ===== TESTING: ======
fprintf('------ TESTING ------ \n')

% Prediction:
[o1,o2,o3] = svmpredict(test_labels, double(test_features), svm_model);

result.model=svm_model;
result.predicted_label=o1;
result.accuracy=o2(1);
result.prob_estimates=o3;
result.options=options;

fprintf('The total classification score is %.4f. \n', ...
    result.accuracy)
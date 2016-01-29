%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% Performs:   CLASSIFICATION                                              %
% on:         MNIST DATASET (can be reduce to some classes only)          %
% using:      SVM.                                                        %
% testing:    Reduced size                                                %
%                                                                         %
% USAGE:                                                                  %
%    1) Change DIR_DATA in the script to point to the dataset             %
%    2) Set the number of training points 'n_training', the number of     %
%    testing points 'n_testing' ('full' for the complete rest of the      %
%    dataset for testing). Set also the scattering architecture (see      %
%    scatnet documentation if needed).                                    %
%    3) Select the label you are insterested in classifying.              %
%                                                                         %
% REQUIRE:                                                                %
%    1) scatnet-master                                                    %
%       http://www.di.ens.fr/data/software/scatnet/                       %
%    2) Libsvm-3.11                                                       %
%       Included in ./Misc just install it                                %
%    3) Download MNIST                                                    %
%       http://yann.lecun.com/exdb/mnist/                                 %
%                                                                         %
% PARAMETERS:                                                             %
%    - Scattering transfrom:                                              %
%       J=3; L=3; M = 2; filter_type='morlet'; oversampling = 2;          %
%    - SVM:                                                               %
%       Cross validation                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


%% ===== Initialization: =====
n_training = 5; % Dataset training + testing has 160 images per class 
n_testing = 20;

rdm_spl = randsample(1:5000, 5000);

rdm_training = rdm_spl(1:n_training);
if strcmp(n_testing,'full')
    rdm_testing = randsample(50000:60000, 10000);
else
    rdm_testing = rdm_spl(n_training+1:n_training+ n_testing);
end

%% ===== PREPARE DATA: =====
% Data path:
im_format = 'ubyte';
dir_data = '/home/jeanbaptiste/Datasets/Mnist/';
addpath(dir_data);

% Loard MNIST:
tmp_mnist = loadMNISTImages('train-images-idx3-ubyte');
tmp_label = loadMNISTLabels('train-labels-idx1-ubyte');

n_label = 10;

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
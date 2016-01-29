%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% Performs:   CLASSIFICATION                                              %
% on:         MNIST DATASET (can be reduce to some classes only)          %
% using:      SCN+SVM.                                                    %
% testing:    Full size                                                   %
%                                                                         %
% USAGE:                                                                  %
%    1) Change DIR_DATA in the script to point to the dataset             %
%    2) Set the number of training points 'n_training', the number of     %
%    testing points 'n_testing' ('full' for the complete rest of the      %
%    dataset for testing). Set also the scattering architecture (see      %
%    scatnet documentation if needed).                                    %
%    3) Select the trained model 'tmp_load' to use for testing            %
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
% Load model
tmp_load = load('./Save/Models/Mnist/SVM/V4_svm_N10_62.mat');
svm_model = tmp_load.svm_model;

% ST Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

%% ===== PREPARE DATA: =====
% Data path:
im_format = 'ubyte';
dir_data = '/home/jeanbaptiste/Datasets/Mnist/';

% Labels available:
S_label = {1, 2 , 3, 4, 5, 6, 7, 8, 9, 0}; 
n_label = length(S_label);

msg_1 = sprintf('Classification between:');
msg_2 = sprintf(' %s',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])

% Generate testing set:  
n_testing = 10000;
rdm_spl = randsample(1:n_testing, n_testing);

%% ===== LOAD TESTSET: =====
addpath(dir_data)
            
% Load MNIST:
tmp_mnist = loadMNISTImages('train-images-idx3-ubyte');
tmp_label = loadMNISTLabels('train-labels-idx1-ubyte');

tmp_mnist= reshape(tmp_mnist, sqrt(size(tmp_mnist,1)), ...
                sqrt(size(tmp_mnist,1)), size(tmp_mnist,2));

test_set = tmp_mnist(:, :, rdm_spl);
test_label = tmp_label(rdm_spl);

%% ===== TESTING: ======
fprintf('------ TESTING ------ \n')

clf_score = 0;
for test_im=1:n_testing
    msg_4 = sprintf('Testing example %i/%i \n', test_im, n_testing);
    msg_5 = sprintf('Clf_score = %i \n',  clf_score);
    fprintf([msg_4, msg_5])
    % Generate testing example:
    x = test_set(:,:,test_im);
    
    Wop = wavelet_factory_2d(size(x), filt_opt, scat_opt);
    S = scat(x, Wop);
    set_S = [{} {S}];
    
    %Reshape ST(set_S{label}{im} into a line vector
    test_feature = [];
    for layer=1:length(S)
        for scale=1:length(S{layer}.signal)
            tmp_vect = reshape(...
                S{layer}.signal{scale}, ...
                1 , numel(S{layer}.signal{scale}));
            test_feature = horzcat(test_feature, tmp_vect);
        end
    end
    % Prediction:
    [o1,o2,o3] = svmpredict(test_label(test_im), double(test_feature), svm_model);
    
    clf_score = clf_score + o2(1)/100;
end
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
% close all


%% ===== Initialization: =====
% Load model
tmp_load = load('./Save/Models/Mnist/IMPORTANT_svm_N5_58.5.mat');
svm_model = tmp_load.svm_model;

% ST Parameters:
filt_opt.J = 4; % scales
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
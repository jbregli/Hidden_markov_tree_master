%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% Performs:   CLASSIFICATION "One vs All"                                 %
% on:         MNIST DATASET (can be reduce to some classes only)          %
% using:      SCN+SVM.                                                    %
% testing:    Reduced size                                                %
%                                                                         %
% USAGE:                                                                  %
%    1) Change DIR_DATA in the script to point to the dataset             %
%    2) Set the number of training points 'n_training', the number of     %
%    testing points 'n_testing' (10000 for the complete testset). Set     %
%    also the scattering architecture (see scatnet documentation if       %
%    needed).                                                             %
%    3) Select the label you are insterested in 'tested_label' (int       %
%    between 0 and 9.                                                     %
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
% 1 vs All testing label:
tested_label = 4;

n_training = 5;
n_testing = 1000; 

% SVM parameters:
svm.model = '-g 0.0018 -c 3';

options.verbose = false;

% ST Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 3; % orientations
scat_opt.M = 2; % layers
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
msg_2 = sprintf(' %i',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])

% Generate training set:        
rdm_spl = randsample(1:5000, 5000);

rdm_training = rdm_spl(1:n_training);

clear allFiles tmp_file tmp_name

%% ===== SCATTERING TRANSFORM: =====
set_S = cell(1,n_label);
for label=1:n_label
    path_to_training = {{dir_data, S_label{label}, im_format}, rdm_training};
    
    fprintf('------ GENERATING CLASS %s %i/%i ------ \n', ...
        S_label{label}, label, n_label)

    % ST:
    set_S{label} = ST_class(path_to_training, filt_opt, scat_opt, n_training, im_format);
end

%% ===== RESHAPE ST:  =====
features = [];
labels = [];
for label=1:n_label
    for im=1:n_training
        %Reshape ST(set_S{label}{im} into a line vector
        tmp_feature = [];
        for layer=1:length(set_S{label}{im})
            for scale=1:length(set_S{label}{im}{layer}.signal)
                tmp_vect = reshape(...
                    set_S{label}{im}{layer}.signal{scale}, ...
                    1 , numel(set_S{label}{im}{layer}.signal{scale}));
                tmp_feature = horzcat(tmp_feature, tmp_vect);
            end
        end
        features = vertcat(features, tmp_feature);
        labels = vertcat(labels, label);
    end
end

% Binarize the problem:
if tested_label == 0
    tested_label = 10;
    labels(labels==0) =10;
end
labels(labels~=tested_label) = 0;
labels(labels~=0) = 1;

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

% Load MNIST:
tmp_mnist = loadMNISTImages('train-images-idx3-ubyte');
tmp_label = loadMNISTLabels('train-labels-idx1-ubyte');

tmp_mnist= reshape(tmp_mnist, sqrt(size(tmp_mnist,1)), ...
                sqrt(size(tmp_mnist,1)), size(tmp_mnist,2));

rdm_spl = randsample(5001:5001 + n_testing, n_testing);
            
test_set = tmp_mnist(:, :, rdm_spl);
test_labels = tmp_label(rdm_spl);

% Binarize the problem:
test_labels(test_labels~=tested_label) = 0;
test_labels(test_labels~=0) = 1;

clf_score = 0;
true_positive = 0;
false_positive = 0;
true_negative = 0;
false_negative = 0;

for test_im=1:n_testing
    msg_4 = sprintf('Testing example %i/%i \n', test_im, n_testing);
    msg_5 = sprintf('Clf_score = %i \n',  clf_score);
    msg_6 = sprintf('TP = %i , FP = %i \n',  true_positive, false_positive);
    msg_7 = sprintf('FN = %i , TN = %i \n',  false_negative, true_negative);

    fprintf([msg_4, msg_5, msg_6, msg_7])
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
    [o1,o2,o3] = svmpredict(test_labels(test_im), double(test_feature), svm_model);

    clf_score = clf_score + o2(1)/100;
    
     if o1 == 1  && o2(1)/100 == 1
         true_positive = true_positive + 1;
     elseif o1 == 1 &&  test_labels(test_im) ~= tested_label
        false_positive = false_positive + 1;
     elseif o1 ~= 1 && o2(1)/100 == 1
        false_negative = false_negative +1;
     elseif o1 ~= 1 && o2(1)/100 ~= 1
        true_negative = true_negative +1;             
     end    
end

fprintf('Clf_score = %i \n',  clf_score);
fprintf('TP = %i , FP = %i \n',  true_positive, false_positive);
fprintf('FN = %i , TN = %i \n',  false_negative, true_negative);

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

options.verbose = false;

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
% aluminium_foil , brown_bread , corduroy , cotton , cracker , linen , 
% orange_peel , sandpaper , sponge , styrofoam
S_label = {1, 2 , 3, 4, 5, 6, 7, 8, 9, 0}; 
n_label = length(S_label);

msg_1 = sprintf('Classification between:');
msg_2 = sprintf(' %s',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])

% Generate training set:        
rdm_spl = randsample(1:5000, 5000);

rdm_training = rdm_spl(1:n_training);
rdm_test = rdm_spl(n_training+1:n_training+ n_testing); % rdm_spl(n_training+1:end);

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

n_test = length(rdm_test);

%% ===== TESTING: ======
fprintf('------ TESTING ------ \n')
n_test = length(rdm_test);

% Generate testing examples: 
for label = 1:n_label
    path_to_test{label} = {{dir_data, S_label{label}, im_format}, ...
        rdm_test};
    S_test{label} = ST_class(path_to_test{label}, ...
        filt_opt, scat_opt, n_test, im_format);
end

test_features = [];
test_labels = [];
for label=1:n_label
    for im=1:n_test
        %Reshape ST(set_S{label}{im} into a line vector
        tmp_feature = [];
        for layer=1:length(S_test{label}{im})
            for scale=1:length(S_test{label}{im}{layer}.signal)
                tmp_vect = reshape(...
                    S_test{label}{im}{layer}.signal{scale}, ...
                    1 , numel(S_test{label}{im}{layer}.signal{scale}));
                tmp_feature = horzcat(tmp_feature, tmp_vect);
            end
        end
        test_features = vertcat(test_features, tmp_feature);
        test_labels = vertcat(test_labels, label);
    end
end

% Prediction:
[o1,o2,o3] = svmpredict(test_labels, double(test_features), svm_model);

result.model=svm_model;
result.predicted_label=o1;
result.accuracy=o2(1);
result.prob_estimates=o3;
result.options=options;

fprintf('The total classification score is %.4f. \n', ...
    result.accuracy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% Performs:   CLASSIFICATION "One vs All"                                 %
% on:         KTH TEXTURE DATASET (can be reduce to some classes only)    %
% using:      SCN+SVM.                                                    %
% testing:    Reduced size                                                %
%                                                                         %
% USAGE:                                                                  %
%    1) Change DIR_DATA in the script to point to the dataset             %
%    2) Set the number of training points 'n_training', the number of     %
%    testing points 'n_testing' ('full' for the complete rest of the      %
%    dataset for testing). Set also the scattering architecture (see      %
%    scatnet documentation if needed).                                    %
%    3) Select the label you are insterested in 'tested_label' (int       %
%    between 1 and 10.                                                    %
%                                                                         %
% REQUIRE:                                                                %
%    1) scatnet-master                                                    %
%       http://www.di.ens.fr/data/software/scatnet/                       %
%    2) Libsvm-3.11                                                       %
%       Included in ./Misc just install it                                %
%    3) Download KTH-TIPS                                                 %
%       http://www.nada.kth.se/cvap/databases/kth-tips/download.html      %
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
tested_label = 4; % 1-10

 % Dataset training + testing has 160 images per class 
n_training = 5;
n_testing = 'full'; 

% ST Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

%% ===== PREPARE DATA: =====
% Data path:
im_format = 'png';
dir_data = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';

S_label = {'aluminium_foil' , 'brown_bread', 'corduroy' , 'cotton' , ...
    'cracker' , 'linen' , 'orange_peel' , 'sandpaper' , 'sponge' ,...
    'styrofoam'}; 
n_label = length(S_label);
tested_label_name = S_label{tested_label};


msg_1 = sprintf('Classification between:');
msg_2 = sprintf(' %s',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])


% Generate training set:
% Create a random sampler:
allFile = zeros(1,n_label);
for label=1:n_label
    tmp_file = dir(fullfile(dir_data, S_label{label}, '/', ['*.' im_format]));
    tmp_name = {tmp_file.name};
    allFile(label) = length(tmp_name);
end
        
n_test = min(allFile);
rdm_spl = randsample(1:n_test, n_test);

rdm_training = rdm_spl(1:n_training);
if strcmp(n_testing, 'full')
    rdm_test = rdm_spl(1:n_test);
else
    rdm_test = rdm_spl(1:n_testing);
    % Training points are reused for testing not perfect:
    % rdm_spl(n_training+1:end);
end

clear allFiles tmp_file tmp_name

%% ===== SCATTERING TRANSFORM: =====
set_S = cell(1,n_label);
for label=1:n_label
    path_to_training = {{dir_data, S_label{label}, im_format}, ...
        rdm_training};
    
    fprintf('------ GENERATING CLASS %s %i/%i ------ \n', ...
        S_label{label}, label, n_label)

    % ST:
    set_S{label} = ST_class(path_to_training, filt_opt, scat_opt, ...
        n_training, im_format);
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

% Binarize the problem:
test_labels(test_labels~=tested_label) = 0;
test_labels(test_labels~=0) = 1;


% Prediction:
[o1,o2,o3] = svmpredict(test_labels, double(test_features), svm_model);

fprintf('One vs All accuracy = %.2f \n', o2(1))
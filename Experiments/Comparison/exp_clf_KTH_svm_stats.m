%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% Performs:   STATISTICS ON CLASSIFICATION                                %
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
%    3) Select the number of experiments you want to do 'n_exp'           %
%    4) Select the name of the csv output 'csv_name'                      %
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

%% ===== Meta parameters: =====
% Number of experiment:
n_exp = 100;
score = cell(1,n_exp);
csv_name = 'stat_clf_KTH_svm_N5_J4_L3.csv';
max_score = 0;

for experiment=1:n_exp
    fprintf('--------- EXPERIMENT %i/%i --------- \n', experiment, n_exp)
    
    %% ===== PREPARE: SCHMT =====
    n_training = 5; % Dataset training + testing has 160 images per class 
    n_testing = 20;

    % ST Parameters:
    filt_opt.J = 4; % scales
    filt_opt.L = 3; % orientations
    filt_opt.filter_type = 'morlet';
    scat_opt.oversampling = 2;
    scat_opt.M = 2;
    
    %% ===== PREPARE DATA: =====
    % Data path:
    im_format = 'png';
    dir_data = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';

    % Labels available:
    % 'aluminium_foil' , 'brown_bread', 'corduroy' , 'cotton' , ...
    %    'cracker' , 'linen' , 'orange_peel' , 'sandpaper' , 'sponge' ,...
    %    'styrofoam'
    S_label = {'aluminium_foil' , 'brown_bread', 'corduroy' , 'cotton' , ...
        'cracker' , 'linen' , 'orange_peel' , 'sandpaper' , 'sponge' ,...
        'styrofoam'}; 
    n_label = length(S_label);

    msg_1 = sprintf('Classification between:');
    msg_2 = sprintf(' %s',  S_label{:});
    msg_3 = sprintf('\n');
    fprintf([msg_1 msg_2 msg_3])


    % Generate training set:
    % Create a random sampler:
    allFile = zeros(1,n_label);
    for label=1:n_label
        tmp_file = dir(fullfile(dir_data, S_label{label}, '/', ...
            ['*.' im_format]));
        tmp_name = {tmp_file.name};
        allFile(label) = length(tmp_name);
    end

    rdm_spl = randsample(1:min(allFile), min(allFile));

    rdm_training = rdm_spl(1:n_training);
    if strcmp(n_testing, 'full')
        rdm_test = rdm_spl(n_training+1:end);
    else
        rdm_test = rdm_spl(n_training+1:n_training+ n_testing);
    end

    clear allFiles tmp_file tmp_name
    
    %% ===== PREPARE: SCORE STRUCTURE =====
    score{experiment}.clf_score = 0.0;
    
    %% ===== SCATTERING TRANSFORM: =====
    set_S = cell(1,n_label);
    for label=1:n_label
        path_to_training = {{dir_data, S_label{label}, im_format}, ...
            rdm_training};

        fprintf('------ GENERATING CLASS %s %i/%i ------ \n', ...
            S_label{label}, label, n_label)

        % ST:
        set_S{label} = ST_class(path_to_training, filt_opt, scat_opt,...
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

    % Prediction:
    [o1,o2,o3] = svmpredict(test_labels, double(test_features), svm_model);

    score{experiment}.clf_score = o2(1);
    
   
    %% ===== CSV WRITTING: =====
    csv_line = [score{experiment}.clf_score];
    
    dlmwrite (csv_name,csv_line, '-append');
    
    clearvarlist = ['clearvarlist'; ...
        setdiff(who,{'n_exp';'score';'csv_name';'max_score'})];
    clear(clearvarlist{:})
end

% Avg clf score:
clf_avg = 0;
for i=1:n_exp
    clf_avg = clf_avg + score{i}.clf_score/n_exp;
end
         
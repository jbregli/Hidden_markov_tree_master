%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% Performs:   STATISTICS ON CLASSIFICATION                                %
% on:         MNIST DATASET (can be reduce to some classes only)          %
% using:      SCN+SVM.                                                    %
% testing:    Reduced size                                                %
%                                                                         %
% USAGE:                                                                  %
%    1) Change DIR_DATA in the script to point to the dataset             %
%    2) Set the number of training points 'n_training', the number of     %
%    testing points 'n_testing'. Set also the scattering architecture     %
%    (see scatnet documentation if needed).                               %
%    3) Select the number of experiments you want to do 'n_exp'           %
%    4) Select the name of the csv output 'csv_name'                      %
%    5) Select the name of the best model name 'mat_name'                 %
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

%% ===== Meta parameters: =====
% Number of experiment:
n_exp = 100;
score = cell(1,n_exp);
csv_name = 'V4_stat_clf_Mnist_svm_N10_J4_L3.csv';
mat_name = 'V4_svm_N10_';
max_score = 0;
format = 'mat';
directory = './Save/Models/Mnist/SVM/';


for experiment=1:n_exp
    fprintf('--------- EXPERIMENT %i/%i --------- \n', experiment, n_exp)
    
    %% ===== PREPARE: SCHMT =====
    n_training = 10; % Dataset training + testing has 160 images per class 
    n_testing = 1000;

    % ST Parameters:
    filt_opt.J = 4; % scales
    filt_opt.L = 3; % orientations
    filt_opt.filter_type = 'morlet';
    scat_opt.oversampling = 6;
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
    rdm_test = rdm_spl(n_training+1:n_training+ n_testing);

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
    
    %% ===== SVM TRAINING: =====
    fprintf('------ TRAINING ------ \n')

    % Shuffle elements
    shuffler = randperm(size(features,1));%
%    
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
    
    fprintf('The total classification score is %.4f. \n', ...
        score{experiment}.clf_score)
    
    if score{experiment}.clf_score < 41
        %delete([directory mat_name '*.mat']);        
        % Save the model:
        save([directory mat_name ...
            num2str(score{experiment}.clf_score) '.mat'],'svm_model');       
    end    
    
    if score{experiment}.clf_score > max_score
        delete([directory mat_name num2str(max_score) '.mat']);        
        max_score = score{experiment}.clf_score;
        save([directory mat_name num2str(max_score) '.mat'],'svm_model');       
    end  

   
    %% ===== CSV WRITTING: =====
    csv_line = [score{experiment}.clf_score];
    
    dlmwrite ([directory csv_name],csv_line, '-append');
    
    clearvarlist = ['clearvarlist'; setdiff(who,{'n_exp';'score';...
        'max_score';'csv_name';'max_score';'format';'directory'; ...
        'mat_name'})];
    clear(clearvarlist{:})
end

% Avg clf score:
clf_avg = 0;
for i=1:n_exp
    clf_avg = clf_avg + score{i}.clf_score/n_exp;
end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test of the STHMT on cropped      % 
% sonar images from MUSSEL AREA C dataset.                                %
% The same experiment is performed n_exp = 30 times to draw statistics    %
% out of the classification score                                         %
% BEST PARAMETER SO FAR: CS=0.8                                           %
% n_image = 0; n_state = 2; n_step = 100; eps_uni= false; cv_sens = 1e-5; %
% filt_opt.J = 5; filt_opt.L = 3; filt_opt.filter_type = 'morlet';        %
% scat_opt.oversampling = 2; scat_opt.M = 2;                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Meta parameter:
% Number of experiment:
n_exp = 100;
score = cell(1,n_exp);
csv_name = 'stat_clf_ripple2_N200_J5_L3.csv';
max_score = 0;
format = 'mat';

for experiment=1:n_exp
    fprintf('--------- EXPERIMENT %i/%i --------- \n', experiment, n_exp)
    
    %% Initialization:
    % Number of image
    n_image = 200; % 0 for all the images

    % Number of states:
    n_state = 2;

    % Number of optimization step:
    n_step = 100;

    % Model distribution:
    distribution = 'MixtGauss';
    % Epsilon uniform over the pixels of a father/son transition
    eps_uni= false;
    % Display error messages:
    verbose = false;
    % Sensibility f the convergence test:
    cv_sens = 1e-5;

    % Data path:
    %dir_training = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Training/';
    dir_training = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops2/Training/';

    % ST Parameters:
    filt_opt.J = 5; % scales
    filt_opt.L = 3; % orientations
    filt_opt.filter_type = 'morlet';
    scat_opt.oversampling = 2;
    scat_opt.M = 2;


    %% Preparation of the structure:
    score{experiment}.clf_score = 0.0;
    score{experiment}.tp_ri_ri = 0.0;
    score{experiment}.fp_ri_se = 0.0;
    score{experiment}.tp_se_se = 0.0;
    score{experiment}.fp_se_ri = 0.0;
    
    %% CLASS 1 - RIPPLE - TRAINING: 
    label_ripple = 'Ripple/'; 
    path_to_training_ripple = fullfile(dir_training, label_ripple);

    fprintf('------ TRAINING RIPPLE ------ \n')

    % ST:
    set_S_ripple = ST_class(path_to_training_ripple, filt_opt, scat_opt, n_image, format);

    % Prepare the scattering structure for HMM:
    for im=1:length(set_S_ripple)
        set_S_ripple{im} = hmm_prepare_S(set_S_ripple{im}, n_state);
    end

    % Hmm model:
    [theta_est_ripple, ~, ~] = ...
        conditional_EM(set_S_ripple, n_step, n_state, distribution, ...
            eps_uni, verbose, 10, cv_sens);

    clear set_S_ripple

    %% CLASS 2 - Mix - TRAINING: 
    label_seabed = 'Seabed'; 
    path_to_training_seabed = fullfile(dir_training, label_seabed);

    fprintf('------ TRAINING SEABED ------ \n')

    % ST:
    set_S_seabed = ST_class(path_to_training_seabed, filt_opt, scat_opt, n_image, format);

    % Prepare the scattering structure for HMM:
    for im=1:length(set_S_seabed)
        set_S_seabed{im} = hmm_prepare_S(set_S_seabed{im}, n_state);
    end

    % Hmm model:
    [theta_est_seabed, ~, ~] = ...
        conditional_EM(set_S_seabed, n_step, n_state, distribution, ...
            eps_uni, verbose, 10, cv_sens);                          

    clear set_S_seabed

    %% MAP - CLASSIFICATION SCORE:
    fprintf('------ TESTING ------ \n')
    n_test = 40;
    TP_ripple = 0; FP_ripple = 0;
    TP_seabed = 0; FP_seabed = 0;

    P_hat_ri_ri = cell(1,n_test);
    P_hat_ri_se = cell(1,n_test);
    P_hat_se_ri = cell(1,n_test);
    P_hat_se_se = cell(1,n_test);

    %dir_test = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Test/';
    dir_test = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops2/Test/';
    
    path_to_test_ripple = fullfile(dir_test, label_ripple);
    S_ripple_test = ST_class(path_to_test_ripple, filt_opt, scat_opt, n_test, format);

    path_to_test_seabed = fullfile(dir_test, label_seabed);
    S_seabed_test = ST_class(path_to_test_seabed, filt_opt, scat_opt, n_test, format);

    % Prepare the scattering structure for HMM:
    for im=1:n_test
        S_ripple_test{im} = hmm_prepare_S(S_ripple_test{im}, n_state);
        S_seabed_test{im} = hmm_prepare_S(S_seabed_test{im}, n_state);

        % MAP model(ripple) = ripple:
        [tmp_P_hat_ri_ri, ~] = ...
            hmm_MAP(S_ripple_test{im}, theta_est_ripple, false);
        [tmp_P_hat_ri_se, ~] = ...
            hmm_MAP(S_ripple_test{im}, theta_est_seabed, false);

        P_hat_ri_ri{im} = mean(mean(tmp_P_hat_ri_ri));
        P_hat_ri_se{im} = mean(mean(tmp_P_hat_ri_se));


        % MAP model(seabed) = seabed:
        [tmp_P_hat_se_ri, ~] = ...
            hmm_MAP(S_seabed_test{im}, theta_est_ripple, false);
        [tmp_P_hat_se_se, ~] = ...
            hmm_MAP(S_seabed_test{im}, theta_est_seabed, false);

        P_hat_se_ri{im} = mean(mean(tmp_P_hat_se_ri));
        P_hat_se_se{im} = mean(mean(tmp_P_hat_se_se));
    end

    for im=1:n_test
        if P_hat_ri_ri{im} > P_hat_ri_se{im}
            TP_ripple = TP_ripple + 1/n_test;
        else
            FP_ripple = FP_ripple + 1/n_test;        
        end

        if P_hat_se_se{im} > P_hat_se_ri{im}
            TP_seabed = TP_seabed + 1/n_test;
        else
            FP_seabed = FP_seabed + 1/n_test;
        end
    end

    score{experiment}.clf_score = (TP_seabed + TP_ripple)/2;
    score{experiment}.tp_ripple = TP_ripple;
    score{experiment}.fp_ripple = FP_ripple;
    score{experiment}.tp_seabed = TP_seabed;
    score{experiment}.fp_seabed = FP_seabed;   
    
    fprintf('The total classification score is %.4f. \n', ...
        score{experiment}.clf_score)
%     fprintf(['TP(ripple = ripple) = %.4f. ', ...
%              'FP(ripple = seabed) = %.4f. \n'], ...
%              score{experiment}.fp_ripple, score{experiment}.fp_ripple)
%     fprintf(['TP(seabed = seabed) = %.4f. ', ...
%              'FP(seabed = ripple) = %.4f. \n'], ...
%              score{experiment}.tp_seabed, score{experiment}.fp_seabed)

    if score{experiment}.clf_score > max_score
        delete('./Saved_model/*.mat');        
        
        max_score = score{experiment}.clf_score;
        % Save the model:
        name_ripple = './Saved_model/theta_est_ripple_';
        name_seabed = './Saved_model/theta_est_seabed_';

        save([name_ripple num2str(max_score) '.mat'],'theta_est_ripple');       
        save([name_seabed num2str(max_score) '.mat'],'theta_est_seabed');  
    end
   
    %% Clearing:
    csv_line = [score{experiment}.clf_score score{experiment}.tp_ripple ...
        score{experiment}.fp_ripple score{experiment}.tp_seabed ...
        score{experiment}.fp_seabed];
    
    dlmwrite (csv_name,csv_line, '-append');
    
    clearvarlist = ['clearvarlist'; setdiff(who,{'n_exp';'score';'csv_name';'max_score';'format'})];
    clear(clearvarlist{:})
end

% Avg clf score:
clf_avg = 0;
for i=1:n_exp
    clf_avg = clf_avg + score{i}.clf_score/n_exp;
end
  
         
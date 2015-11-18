%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test on few classes from the KTH  %
% texture dataset.                                                        %
%                                                                         %
% BEST PARAMETER SO FAR: CS=0.8                                           %
% n_image = 0; n_state = 2; n_step = 100; eps_uni= false; cv_sens = 1e-5; %
% filt_opt.J = 5; filt_opt.L = 3; filt_opt.filter_type = 'morlet';        %
% scat_opt.oversampling = 2; scat_opt.M = 2;                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% ===== Initialization: =====
n_training = 60; % Dataset training + testing has 81 images per class 

% EM parameters:
EM_metaparameters.n_step = inf;
EM_metaparameters.n_state = 2;
EM_metaparameters.distribution = 'MixtGauss';
EM_metaparameters.eps_uni = true;
EM_metaparameters.mixing = 10;
EM_metaparameters.cv_sens = 1e-5;
EM_metaparameters.cv_steps = 7;
EM_metaparameters.cv_ratio = 0.98;

options.verbose = false;

% ST Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 4; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

%% ===== PREPARE DATA: =====
% Data path:
im_format = 'png';
dir_data = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';

% Labels available:
% aluminium_foil , brown_bread , corduroy , cotton , cracker , linen , 
% orange_peel , sandpaper , sponge , styrofoam
label_1 = 'corduroy';
label_2 = 'brown_bread'; 

msg_1 = sprintf('Classification between %s and %s. \n',  label_1, label_2);
fprintf(msg_1)


% Generate training set:
% Create a random sampler:
allFiles_1 = dir(fullfile(dir_data, label_1, '/', ['*.' im_format]));
allNames_1 = {allFiles_1.name};
allFiles_2 = dir(fullfile(dir_data, label_2, '/', ['*.' im_format]));
allNames_2 = {allFiles_2.name};
        
rdm_spl = randsample(1:min(length(allNames_1),length(allNames_2)), ...
    min(length(allNames_1),length(allNames_2)));

rdm_training = rdm_spl(1:n_training);
rdm_test = rdm_spl(n_training+1:end);

clear allFiles_1 allNames_1 allFiles_2 allNames_2

%% ===== CLASS 1 - TRAINING: =====
path_to_training_c1 = {{dir_data, label_1, im_format}, rdm_training};

fprintf('------ TRAINING CLASS 1 ------ \n')

% ST:
set_S_c1 = ST_class(path_to_training_c1, filt_opt, scat_opt, n_training, im_format);

% Prepare the scattering structure for HMM:
for im=1:length(set_S_c1)
    set_S_c1{im} = hmm_prepare_S(set_S_c1{im}, EM_metaparameters.n_state);
end

% Hmm model:
[theta_est_c1, cv_stat_c1, dob1] = ...
    conditional_EM(set_S_c1, EM_metaparameters, options);

clear set_S_c1

%% CLASS 2 - TRAINING: 
path_to_training_c2 = {{dir_data, label_2, im_format}, rdm_training};

fprintf('------ TRAINING CLASS 2 ------ \n')

% ST:
set_S_c2 = ST_class(path_to_training_c2, filt_opt, scat_opt, n_training, im_format);

% Prepare the scattering structure for HMM:
for im=1:length(set_S_c2)
    set_S_c2{im} = hmm_prepare_S(set_S_c2{im}, EM_metaparameters.n_state);
end

% Hmm model:
[theta_est_c2,  cv_stat_c2, dob2] = ...
    conditional_EM(set_S_c2, EM_metaparameters, options);                          

clear set_S_c2
    
%% MAP - CLASSIFICATION SCORE:
fprintf('------ TESTING ------ \n')
n_test = length(rdm_test);
TP_c1 = 0; FP_c1 = 0;
TP_c2 = 0; FP_c2 = 0;

P_hat_c1_c1 = cell(1,n_test);
P_hat_c1_c2 = cell(1,n_test);
P_hat_c2_c1 = cell(1,n_test);
P_hat_c2_c2 = cell(1,n_test);


path_to_test_c1 = {{dir_data, label_1, im_format}, rdm_test};
S_test_c1 = ST_class(path_to_test_c1, filt_opt, scat_opt, n_test, im_format);

path_to_test_c2 = {{dir_data, label_2, im_format}, rdm_test};
%S_seabed_test = ST_class(path_to_test_seabed, filt_opt, scat_opt, n_test, format);
S_test_c2 = ST_class(path_to_test_c2, filt_opt, scat_opt, n_test, im_format);


% Prepare the scattering structure for HMM:
for im=1:n_test
    S_test_c1{im} = hmm_prepare_S(S_test_c1{im}, EM_metaparameters.n_state);
    S_test_c2{im} = hmm_prepare_S(S_test_c2{im}, EM_metaparameters.n_state);
    
    % MAP model(ripple) = ripple:
    [tmp_P_hat_c1_c1, ~] = ...
        hmm_MAP(S_test_c1{im}, theta_est_c1, false);
    [tmp_P_hat_c1_c2, ~] = ...
        hmm_MAP(S_test_c1{im}, theta_est_c2, false);
    
    P_hat_c1_c1{im} = mean(mean(tmp_P_hat_c1_c1));
    P_hat_c1_c2{im} = mean(mean(tmp_P_hat_c1_c2));

    
    % MAP model(seabed) = seabed:
    [tmp_P_hat_c2_c1, ~] = ...
        hmm_MAP(S_test_c2{im}, theta_est_c1, false);
    [tmp_P_hat_c2_c2, ~] = ...
        hmm_MAP(S_test_c2{im}, theta_est_c2, false);
    
    P_hat_c2_c1{im} = mean(mean(tmp_P_hat_c2_c1));
    P_hat_c2_c2{im} = mean(mean(tmp_P_hat_c2_c2));

end

for im=1:n_test  
    if P_hat_c1_c1{im} > P_hat_c1_c2{im}
        TP_c1 = TP_c1 + 1/n_test;
    else
        FP_c1 = FP_c1 + 1/n_test;        
    end
    
    if P_hat_c2_c2{im} > P_hat_c2_c1{im}
        TP_c2 = TP_c2 + 1/n_test;
    else
        FP_c2 = FP_c2 + 1/n_test;
    end
end

fprintf('The total classification score is %.4f. \n', ...
    (TP_c2 + TP_c1)/2)
fprintf(['TP(%s = %s) = %.4f. ', ...
         'FP(%s = %s) = %.4f. \n'], ...
         label_1, label_1, TP_c1, label_1, label_2, FP_c1)
fprintf(['TP(%s = %s) = %.4f. ', ...
         'FP(%s = %s) = %.4f. \n'], ...
         label_2, label_2, TP_c2, label_2, label_1, FP_c2)
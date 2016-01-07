%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test of the STHMT on genereted    %
% cropped sonar images from MUSSEL AREA C dataset.                        %
%                                                                         %
% BEST PARAMETER SO FAR: CS=0.8                                           %
% n_image = 0; n_state = 2; n_step = 100; eps_uni= false; cv_sens = 1e-5; %
% filt_opt.J = 5; filt_opt.L = 3; filt_opt.filter_type = 'morlet';        %
% scat_opt.oversampling = 2; scat_opt.M = 2;                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% ===== Initialization: =====
n_image = 100; % 0 for all the images

% EM parameters:
EM_metaparameters.n_step = 40;
EM_metaparameters.n_state = 2;
EM_metaparameters.distribution = 'MixtGauss';
EM_metaparameters.eps_uni = true;
EM_metaparameters.mixing = 10;
EM_metaparameters.cv_sens = 1e-4;
EM_metaparameters.cv_steps = 7;
EM_metaparameters.cv_ratio = 0.98;

options.verbose = false;

% ST Parameters:
filt_opt.J = 5; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;


%% ===== PREPARE DATA: =====
% Data path:
im_format = 'png';
dir_training = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Training/';

label_ripple = 'Ripple/'; 
label_seabed = 'Seabed'; 

%% ===== CLASS 1 - RIPPLE - TRAINING: =====
path_to_training_ripple = fullfile(dir_training, label_ripple);

fprintf('------ TRAINING RIPPLE ------ \n')

% ST:
set_S_ripple = ST_class(path_to_training_ripple, filt_opt, scat_opt, n_image, im_format);

% Prepare the scattering structure for HMM:
for im=1:length(set_S_ripple)
    set_S_ripple{im} = hmm_prepare_S(set_S_ripple{im}, EM_metaparameters.n_state);
end

% Hmm model:
[theta_est_ripple, ~, ~] = ...
    conditional_EM(set_S_ripple, EM_metaparameters, options);

clear set_S_ripple

%% ===== CLASS 2 - Mix - TRAINING: =====
path_to_training_seabed = fullfile(dir_training, label_seabed);

fprintf('------ TRAINING SEABED ------ \n')

% ST:
set_S_seabed = ST_class(path_to_training_seabed, filt_opt, scat_opt, n_image, im_format);

% Prepare the scattering structure for HMM:
for im=1:length(set_S_seabed)
    set_S_seabed{im} = hmm_prepare_S(set_S_seabed{im}, EM_metaparameters.n_state);
end

% Hmm model:
[theta_est_seabed, ~, ~] = ...
    conditional_EM(set_S_seabed, EM_metaparameters, options);                          

clear set_S_seabed
    
%% ===== MAP - CLASSIFICATION SCORE: =====
fprintf('------ TESTING ------ \n')
n_test = 40;
TP_ripple = 0; FP_ripple = 0;
TP_seabed = 0; FP_seabed = 0;

P_hat_ri_ri = cell(1,n_test);
P_hat_ri_se = cell(1,n_test);
P_hat_se_ri = cell(1,n_test);
P_hat_se_se = cell(1,n_test);

%dir_test = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Test/';
dir_test = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Test/';

path_to_test_ripple = fullfile(dir_test, label_ripple);
S_ripple_test = ST_class(path_to_test_ripple, filt_opt, scat_opt, n_test, im_format);

path_to_test_seabed = fullfile(dir_test, label_seabed);
%S_seabed_test = ST_class(path_to_test_seabed, filt_opt, scat_opt, n_test, im_format);
S_seabed_test = ST_class(path_to_test_seabed, filt_opt, scat_opt, n_test, im_format);

% +++ Mean P for normalization:
sum_P_ri = 0;
sum_P_se = 0;


% Prepare the scattering structure for HMM:
for im=1:n_test
    S_ripple_test{im} = hmm_prepare_S(S_ripple_test{im}, EM_metaparameters.n_state);
    S_seabed_test{im} = hmm_prepare_S(S_seabed_test{im}, EM_metaparameters.n_state);
    
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
    
    P_hat_se_ri{im} = mean(mean(exp(tmp_P_hat_se_ri)));
    P_hat_se_se{im} = mean(mean(exp(tmp_P_hat_se_se)));

%     % +++
%     sum_P_ri = sum_P_ri + P_hat_ri_ri{im};
%     sum_P_se = sum_P_se + P_hat_se_se{im};
end

for im=1:n_test
    % +++
%     P_hat_ri_ri{im} = (P_hat_ri_ri{im} * n_test) / sum_P_ri;
%     P_hat_ri_se{im} = (P_hat_ri_se{im} * n_test) / sum_P_se;
%     
%     P_hat_se_se{im} = (P_hat_se_se{im} * n_test) / sum_P_se;
%     P_hat_se_ri{im} = (P_hat_se_ri{im} * n_test) / sum_P_ri;    
    %+++
    
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

fprintf('The total classification score is %.4f. \n', ...
    (TP_seabed + TP_ripple)/2)
fprintf(['TP(ripple = ripple) = %.4f. ', ...
         'FP(ripple = seabed) = %.4f. \n'], ...
         TP_ripple, FP_ripple)
fprintf(['TP(seabed = seabed) = %.4f. ', ...
         'FP(seabed = ripple) = %.4f. \n'], ...
         TP_seabed, FP_seabed)
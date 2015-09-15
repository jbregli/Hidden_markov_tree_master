%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test of the STHMT on genereted    %
% cropped sonar images from MUSSEL AREA C dataset.                        %
% The SCHMT is trained to model squares.                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Initialization:
% Number of states:
n_state = 2;

% Number of optimization step:
n_step = 100;

% Model distribution:
distribution = 'MixtGauss';
% Epsilon uniform over the pixels of a father/son transition
eps_uni= true;
% Display error messages:
verbose = false;
% Sensibility f the convergence test:
cv_sens = 1e-6;

% Data path:
dir_training = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Training/';

% ST Parameters:
filt_opt.J = 4; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

%% CLASS 1 - RIPPLE - TRAINING: 
label = 'Ripple/'; 
path_to_training_ripple = fullfile(dir_training, label);

fprintf('------ TRAINING RIPPLE ------ \n')

% ST:
set_S_ripple = ST_class(path_to_training_ripple, filt_opt, scat_opt);

% Prepare the scattering structure for HMM:
for im=1:length(set_S_ripple)
    set_S_ripple{im} = hmm_prepare_S(set_S_ripple{im}, n_state);
end

% Hmm model:
[theta_est_ripple, ~, ~] = ...
    conditional_EM(set_S_ripple, n_step, n_state, distribution, ...
        eps_uni, verbose, 10, cv_sens);


%% CLASS 1 - CIRCLE - TRAINING: 
label = 'Mix'; 
path_to_training_seabed = fullfile(dir_training, label);

fprintf('------ TRAINING SEABED ------ \n')

% ST:
set_S_seabed = ST_class(path_to_training_seabed, filt_opt, scat_opt);

% Prepare the scattering structure for HMM:
for im=1:length(set_S_seabed)
    set_S_seabed{im} = hmm_prepare_S(set_S_seabed{im}, n_state);
end

% Hmm model:
[theta_est_seabed, ~, ~] = ...
    conditional_EM(set_S_seabed, n_step, n_state, distribution, ...
        eps_uni, verbose, 10, cv_sens);                          

%% MAP - CLASSIFICATION SCORE:
fprintf('------ TESTING ------ \n')
n_test = 40;
score_ripple = 0;
score_seabed = 0;

P_hat_ri_ri = cell(1,n_test);
P_hat_ri_se = cell(1,n_test);
P_hat_se_ri = cell(1,n_test);
P_hat_se_se = cell(1,n_test);

dir_test = '/home/jeanbaptiste/Datasets/Sonar/Area_C_crops/Test/';

label = 'Ripple/'; 
path_to_test_ripple = fullfile(dir_test, label);
S_ripple_test = ST_class(path_to_test_ripple, filt_opt, scat_opt);

label = 'Mix/'; 
path_to_test_seabed = fullfile(dir_test, label);
S_seabed_test = ST_class(path_to_test_seabed, filt_opt, scat_opt);

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

    if max(P_hat_ri_se{im}, P_hat_ri_ri{im}) == P_hat_ri_ri{im}
        score_ripple = score_ripple + 1/n_test;
    end
    
    % MAP model(seabed) = seabed:
    [tmp_P_hat_se_ri, ~] = ...
        hmm_MAP(S_seabed_test{im}, theta_est_ripple, false);
    [tmp_P_hat_se_se, ~] = ...
        hmm_MAP(S_seabed_test{im}, theta_est_seabed, false);
    
    P_hat_se_ri{im} = mean(mean(tmp_P_hat_se_ri));
    P_hat_se_se{im} = mean(mean(tmp_P_hat_se_se));

    if max(P_hat_se_se{im}, P_hat_se_ri{im}) == P_hat_se_se{im}
        score_seabed = score_seabed + 1/n_test;
    end
end

fprintf('The total classification score is %.4f. \n', ...
    (score_seabed + score_ripple)/2)
fprintf('The classification score for "ripple" is %.4f. \n', ...
    score_ripple)
fprintf('The total classification score for "seabed" is %.4f. \n', ...
    score_seabed)
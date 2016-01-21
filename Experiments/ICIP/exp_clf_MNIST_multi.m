%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a classification test on few classes from the      %
% MNIST dataset.                                                          %
%                                                                         %
% BEST PARAMETER SO FAR: CS=0.6                                           %
% --- ST ---                                                              %                                                   
% filt_opt.J = 4; filt_opt.L = 3; filt_opt.filter_type = 'morlet';        %
% scat_opt.oversampling = 2; scat_opt.M = 2;                              %
% --- EM ---                                                              %                                                   
% EM_meta.n_step = 50; EM_meta.n_state = 2;                               %
% EM_meta.distribution = 'MixtGauss'; EM_meta.eps_uni = true;             %
% EM_meta.mixing = 10; EM_meta.cv_sens = 1e-4; EM_meta.cv_steps = 5;      %
% EM_meta.cv_ratio = 0.8; EM_meta.rerun = true; EM_meta.rerun_count = 0;  %
% EM_meta.rerun_lim = 20;                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% close all

%% ===== Initialization: =====
n_training = 5; % Dataset training + testing has 160 images per class 

% EM parameters:
EM_meta.n_step = 50;
EM_meta.n_state = 2;
EM_meta.distribution = 'MixtGauss';
EM_meta.eps_uni = false;
EM_meta.mixing = 10;
EM_meta.cv_sens = 1e-4;
EM_meta.cv_steps = 5;
EM_meta.cv_ratio = 0.8;
EM_meta.rerun = true;
EM_meta.rerun_count = 0;
EM_meta.rerun_lim = 20;

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
S_label = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0}; 
n_label = length(S_label);

msg_1 = sprintf('Classification between:');
msg_2 = sprintf(' %s',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])

% Generate training set:        
rdm_spl = randsample(1:5000, 5000);

rdm_training = rdm_spl(1:n_training);
rdm_test = rdm_spl(n_training+1:n_training+ 20); % rdm_spl(n_training+1:end);

clear allFiles tmp_file tmp_name

%% ===== TRAINING: =====
theta_est = cell(1,n_label);
cv_stat = cell(1,n_label);
for label=1:n_label
    path_to_training = {{dir_data, S_label{label}, im_format}, rdm_training};
    
    fprintf('------ TRAINING CLASS %s %i/%i ------ \n', ...
        S_label{label}, label, n_label)

    % ST:
    set_S = ST_class(path_to_training, filt_opt, scat_opt, n_training, im_format);

    % Prepare the scattering structure for HMM:
    for im=1:length(set_S)
        set_S{im} = hmm_prepare_S(set_S{im}, EM_meta.n_state);
    end

    % Hmm model:
    [theta_est{label}, cv_stat{label}, ~] = ...
        conditional_EM(set_S, EM_meta, options);
end
% tmp_load = load('./Save/Models/Mnist/IMPORTANT_theta_N5_est_0.8.mat');
% theta_est = tmp_load.theta_est;


%% MAP - CLASSIFICATION SCORE:
fprintf('------ TESTING ------ \n')

n_test = length(rdm_test);

% confusion = matrix(n_label, n_label)
% confusion(i,j): number of object i labeled as j
confusion = zeros(n_label, n_label);

% P_hat = matrix(n_label, n_label)
% P_hat(i,j) probability of an object of class i given model for class j
P_hat = zeros(n_test * n_label, n_label);
test_label = zeros(n_test * n_label, 1);

path_to_test = cell(1, n_label);
S_test = cell(1, n_label);

for label = 1:n_label
    path_to_test{label} = {{dir_data, S_label{label}, im_format}, ...
        rdm_test};
    S_test{label} = ST_class(path_to_test{label}, ...
        filt_opt, scat_opt, n_test, im_format);
end

% Prepare the scattering structure for HMM:
im_count = 1;
for im=1:n_test
    for label=1:n_label
        S_test{label}{im} = hmm_prepare_S(S_test{label}{im},...
            EM_meta.n_state);
        
        for label_model=1:n_label
            [tmp_P_hat, ~] = ...
                hmm_MAP(S_test{label}{im}, theta_est{label_model}, false);
            
            P_hat(im_count, label_model) = mean(mean(tmp_P_hat));
            test_label(im_count) = label;
        end
        im_count = im_count + 1;

        % Rescaling:
%         tmp_mean = mean(P_hat(:,:, im_count));
%         P_hat(:,:, im_count) = P_hat(:,:, im_count) ./ repmat(tmp_mean,length(tmp_mean),1);         
    end
end
for im=1:(n_test * n_label)
    [~, max_idx] = max(P_hat(im,:));
    if max_idx == test_label(im)
        confusion(test_label(im),test_label(im)) = ...
            confusion(test_label(im),test_label(im)) + 1/n_test;
    else
    	confusion(test_label(im), max_idx) = ...
            confusion(test_label(im), max_idx) + 1/n_test;
	end
end

fprintf('The total classification score is %.4f. \n', ...
    sum(diag(confusion)) / n_label)
% fprintf(['TP(%s = %s) = %.4f. ', ...
%          'FP(%s = %s) = %.4f. \n'], ...
%          label_1, label_1, TP_c1, label_1, label_2, FP_c1)
% fprintf(['TP(%s = %s) = %.4f. ', ...
%          'FP(%s = %s) = %.4f. \n'], ...
%          label_2, label_2, TP_c2, label_2, label_1, FP_c2)
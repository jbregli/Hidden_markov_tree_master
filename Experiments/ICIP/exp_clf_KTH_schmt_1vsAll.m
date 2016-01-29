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
% close all


%% ===== Initialization: =====
n_testing = 810;
% 1 vs All testing label:
tested_label = 'styrofoam';
% Labels available:
% 'aluminium_foil' , 'brown_bread', 'corduroy' , 'cotton' , ...
%    'cracker' , 'linen' , 'orange_peel' , 'sandpaper' , 'sponge' ,...
%    'styrofoam'

% Load model
tmp_load = load('./Save/Models/KTH/V4_theta_N2_est_0.1654.mat');
theta_est = tmp_load.theta_est;

% ST Parameters:
filt_opt.J = 3; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

% EM parameters:
EM_meta.n_step = 50;
EM_meta.n_state = 2;
EM_meta.distribution = 'MixtGauss';
EM_meta.eps_uni = true;
EM_meta.mixing = 10;
EM_meta.cv_sens = 1e-4;
EM_meta.cv_steps = 5;
EM_meta.cv_ratio = 0.8;
EM_meta.rerun = true;
EM_meta.rerun_count = 0;
EM_meta.rerun_lim = 20;

%% ===== PREPARE DATA: =====
% Data path:
im_format = 'png';
dir_data = '/home/jeanbaptiste/Datasets/Texture/KTH_TIPS/';

S_label = {'aluminium_foil' , 'brown_bread', 'corduroy' , 'cotton' , ...
        'cracker' , 'linen' , 'orange_peel' , 'sandpaper' , 'sponge' ,...
        'styrofoam'}; 
n_label = length(S_label);

msg_1 = sprintf('Classification between:');
msg_2 = sprintf(' %s',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])


% Generate file index:
test_data = {};
test_label = {};
allFile = zeros(1,n_label);
im_count = 1;
for label=1:n_label
    tmp_file = dir(fullfile(dir_data, S_label{label}, '/', ['*.' im_format]));
    tmp_name = {tmp_file.name};
    allFile(label) = length(tmp_name);  
    for i=1:allFile(label)
        test_data{im_count} = ...
            fullfile(dir_data, S_label{label}, '/', tmp_name(i));
        test_label{im_count} = S_label{label};
        im_count = im_count + 1;
    end
end

n_test = min(allFile);
rdm_test = randsample(1:n_test, n_test);

% clear allFiles tmp_file tmp_name



%% ===== TESTING: ======
fprintf('------ TESTING ------ \n')

clf_score = 0;
confusion = zeros(n_label);
for test_im=1:n_testing
    msg_4 = sprintf('Testing example %i/%i \n', test_im, n_testing);
    msg_5 = sprintf('Clf_score = %i \n',  clf_score);
    
    if test_im == 1
        fprintf([msg_4, msg_5])
        reverseStr = repmat(sprintf('\b'), 1, ....
            length(msg_4) + length(msg_5));
    else
        fprintf([reverseStr, msg_4, msg_5])
        reverseStr = repmat(sprintf('\b'), 1, ....
            length(msg_4) + length(msg_5));
    end
    
    % Generate testing example: 
	x = im2double(imread(char(test_data{test_im})));
    % +++ Resize to avoid wrong image sizes:
    x = x(50:100,50:100);
    
    Wop = wavelet_factory_2d(size(x), filt_opt, scat_opt);
    S = scat(x, Wop);

    S_test = hmm_prepare_S(S, EM_meta.n_state);
            
    % Prediction:
    P_hat = zeros(n_label,1);
	for label_model=1:n_label
    	[tmp_P_hat, ~] = ...
        	hmm_MAP(S_test, theta_est{label_model}, false);
        P_hat(label_model) = mean(mean(tmp_P_hat));
	end
        
    % Rescaling:
    tmp_mean = mean(P_hat);
    P_hat = P_hat ./ tmp_mean;
     
    % Classification
    [~,max_idx] = max(P_hat);
    
    if strcmp(S_label{max_idx},tested_label) ...
            && strcmp(test_label(test_im),tested_label)
        clf_score = clf_score + 1;
    elseif not(strcmp(S_label{max_idx},tested_label)) ...
            &&  not(strcmp(test_label(test_im),tested_label))
        clf_score = clf_score + 1;
    end   
end
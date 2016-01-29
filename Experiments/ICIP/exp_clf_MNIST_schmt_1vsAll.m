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
% Load model
tmp_load = load('./Save/Models/Mnist/V4_theta_N5_est_0.5.mat');
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
im_format = 'ubyte';
dir_data = '/home/jeanbaptiste/Datasets/Mnist/';

% Labels available:
S_label = {1, 2 , 3, 4, 5, 6, 7, 8, 9, 0}; 
n_label = length(S_label);

msg_1 = sprintf('Classification between:');
msg_2 = sprintf(' %s',  S_label{:});
msg_3 = sprintf('\n');
fprintf([msg_1 msg_2 msg_3])

% Generate testing set:  
n_testing = 1000;
rdm_spl = randsample(1:n_testing, n_testing);

% 1 vs All testing label:
tested_label = 0;

%% ===== LOAD TESTSET: =====
addpath(dir_data)
            
% Load MNIST:
tmp_mnist = loadMNISTImages('train-images-idx3-ubyte');
tmp_label = loadMNISTLabels('train-labels-idx1-ubyte');

tmp_mnist= reshape(tmp_mnist, sqrt(size(tmp_mnist,1)), ...
                sqrt(size(tmp_mnist,1)), size(tmp_mnist,2));

test_set = tmp_mnist(:, :, rdm_spl);
test_label = tmp_label(rdm_spl);

%% ===== TESTING: ======
fprintf('------ TESTING ------ \n')

true_positive = 0;
false_positive = 0;
true_negative = 0;
false_negative = 0;

confusion = zeros(n_label);
for test_im=1:n_testing
    msg_4 = sprintf('Testing example %i/%i \n', test_im, n_testing);
    msg_5 = sprintf('TP = %i , FP = %i \n',  true_positive, false_positive);
    msg_6 = sprintf('FN = %i , TN = %i \n',  false_negative, true_negative);
   
    if test_im ==1
        fprintf([msg_4, msg_5, msg_6])
    else
        fprintf([reverseStr, msg_4, msg_5, msg_6])
    end
    reverseStr = repmat(sprintf('\b'), 1, ....
    	length(msg_4) + length(msg_5) + length(msg_6));
    % Generate testing example: 
    x = test_set(:,:,test_im);
    
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
    
    if tested_label == 0
        if max_idx == 10 &&  test_label(test_im) == 0
            true_positive = true_positive + 1;
        elseif max_idx == 10 &&  test_label(test_im) ~= 0
            false_positive = false_positive + 1;
        elseif max_idx ~= 10 && test_label(test_im) == 0
            false_negative = false_negative +1;
        elseif max_idx ~= 10 && test_label(test_im) ~= 0
            true_negative = true_negative +1;
        end
    else
        if max_idx == tested_label && test_label(test_im) == tested_label
            true_positive = true_positive + 1;
        elseif max_idx == tested_label &&  test_label(test_im) ~= tested_label
            false_positive = false_positive + 1;
        elseif max_idx ~= tested_label && test_label(test_im) == tested_label
            false_negative = false_negative +1;
        elseif max_idx ~= tested_label && test_label(test_im) ~= tested_label
            true_negative = true_negative +1;             
        end
    end   
end

fprintf('TP = %i , FP = %i \n',  true_positive, false_positive)
fprintf('FN = %i , TN = %i \n',  false_negative, true_negative);
fprintf('TOT true= %i , TOT false = %i \n', ...
    sum(test_label==tested_label), sum(test_label~=tested_label))

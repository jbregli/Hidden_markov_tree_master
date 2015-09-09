%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script first realize the scattering transform of an images.        %
% Then it computes the correlation between each father node and its       %
% children and display the interesting information.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:

%% Image:
% Easy:
% x = generate_1st_diag(false);
% x = generate_square(false, false); % empty, noise
% temp = im2double(rgb2gray(...
%               imread('/home/jeanbaptiste/Datasets/Sonar/UDRC_datacentre_MCM_sonar_data/Area_B/MUSCLE_COL2_080429_2_7_s_6804_6969_40_150.png')));
% x = temp(1168:1808,1677:2317);

% Medium:
% x = generate_circle(false, true); % empty, noise

% Hard:
x = uiuc_sample;
% x = mandrill;
% x = lena;



%% Scattering transform:
% Parameters:
filt_opt.J = 4; % scales
filt_opt.L = 4; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 3;
% filt_opt = struct();
% scat_opt = struct();


% LOOP OVER THE IMAGES:
% Precompute the WT op thatwill be applied to the image:
[Wop, filters] = wavelet_factory_2d(size(x),filt_opt, scat_opt);

% ST:
[S, U] = scat(x, Wop);

%% Display:
%image_scat(S, true, true);

%% Correlation:
corr_strct = correlation_fc_allNet(S, true, false); % plot, verbose
mean_corr = mean([corr_strct.corr{:}]);
std_corr = std([corr_strct.corr{:}]);

disp(['Average correlation between a father and a child = ' num2str(mean_corr)])
disp(['Std correlation between a father and a child = ' num2str(std_corr)])
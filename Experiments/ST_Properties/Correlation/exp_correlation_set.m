%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script first realize the scattering transform of several images.   %
% Then it computes the correlation between each father node and its       %
% children and display the interesting information.                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:
testset = {};

%% Image:
% Easy:
testset.diag = {};
testset.diag.x = generate_1st_diag(false);

testset.square = {};
testset.square.x = generate_square(false, false); % empty, noise

% testset.udrcB = {};
% temp = im2double(rgb2gray(...
%               imread('/home/jeanbaptiste/Datasets/Sonar/UDRC_datacentre_MCM_sonar_data/Area_B/MUSCLE_COL2_080429_2_7_s_6804_6969_40_150.png')));
% testset.udrcB.x = temp(1168:1808,1677:2317);

% Medium:
testset.circle = {};
testset.circle.x = generate_circle(false, true); % empty, noise

% Hard:
testset.uiuc = {};
testset.uiuc.x = uiuc_sample;

testset.mandrill = {};
testset.mandrill.x = mandrill;

testset.lena = {};
testset.lena.x = lena;


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
fields = fieldnames(testset);
for i = 1:numel(fields)
    disp(['--- ' fields{i} ' ---'])
    % Precompute the WT op thatwill be applied to the image:
    [Wop, filters] = wavelet_factory_2d(size(testset.(fields{i}).x),...
                                                       filt_opt, scat_opt);
                                                   
    % ST:
    [S, U] = scat(testset.(fields{i}).x, Wop);

    %% Display:
    %image_scat(S, true, true);

    %% Correlation:
    corr_strct = correlation_fc_allNet(S);
    testset.(fields{i}).mean_corr = mean([corr_strct.corr{:}]);
    testset.(fields{i}).std_corr = std([corr_strct.corr{:}]);

    disp(['    Average correlation between a father and a child = ' ...
                                   num2str(testset.(fields{i}).mean_corr)])
    disp(['    Std of correlation between a father and a child = ' ...
                                   num2str(testset.(fields{i}).std_corr)])
end
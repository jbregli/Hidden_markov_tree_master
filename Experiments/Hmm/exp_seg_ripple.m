%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script realizes a segmentation test of the STHMT on genereted      %
% cropped sonar images from MUSSEL AREA C dataset.                        %
% The SCHMT is trained to model squares.                                  %
%                                                                         %
% BEST PARAMETER SO FAR: CS=0.8                                           %
% n_image = 0; n_state = 2; n_step = 100; eps_uni= false; cv_sens = 1e-5; %
% filt_opt.J = 5; filt_opt.L = 3; filt_opt.filter_type = 'morlet';        %
% scat_opt.oversampling = 2; scat_opt.M = 2;                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% === load data -------------------------------------------------------
dname = '/home/jeanbaptiste/Datasets/Sonar/UDRC_datacentre_MCM_sonar_data/Area_C/';

files = dir( fullfile(dname,'*.png'));
flist = {files.name};
n_file =  length(flist);
%file_ite = randperm(n_file);
%fname = flist{file_ite(1)};
 fname = flist{40}; %---> ok results
dfname = strcat(dname, fname);

Im = im2double(rgb2gray(imread(dfname)));
Im = Im(850:end -850, 900:end - 900);

%imagesc(Im)
colormap gray

imwrite(Im,'./Save/test_segmentation.png','png')

%%

% Number of states:
n_state = 2;

% Model distribution:
distribution = 'MixtGauss';
% Epsilon uniform over the pixels of a father/son transition
eps_uni= false;
% Display error messages:
verbose = false;
% Sensibility f the convergence test:
cv_sens = 1e-5;


% ST Parameters:
filt_opt.J = 5; % scales
filt_opt.L = 3; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

% Load model:
theta_est_ripple = load('./Save/Saved_model/theta_est_ripple_0.9.mat');
theta_est_seabed = load('./Save/Saved_model/theta_est_seabed_0.9.mat');

theta_est_ripple = theta_est_ripple.theta_est_ripple;
theta_est_seabed = theta_est_seabed.theta_est_seabed;

%% MAP - SEGMENTATION SCORE:
fprintf('------ TESTING ------ \n')
dname = './Save/';
fname_seg = 'test_segmentation.png';
dfname = strcat(dname, fname_seg);
X = im2double(imread(dfname));

s_X = size(X);
s_X(:) = floor(s_X(:)/100) * 100 + 1 ;

x_range = 1:100:s_X(1);
y_range = 1:100:s_X(2);

segmentation = zeros(s_X(1),s_X(2));
pmap_ripple = zeros(s_X(1),s_X(2));
pmap_seabed = zeros(s_X(1),s_X(2));

reverseStr = '';
counter= 1;
n_patch = (length(y_range)-1)*(length(x_range)-1);



for i=1:(length(x_range)-1)
    for j=1:(length(y_range)-1)
        
        % Print time remaining:
        if counter == 1
            tic
            msg = sprintf('--- Patch %i/%i', counter, ...
                (length(y_range)-1)*(length(x_range)-1));
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));   
        elseif counter == 2
           time = toc;
           
           msg = sprintf(['--- Patch %i/%i --- Expected remaining ' ...
                'time: %.4f s. \r '], ...
                counter, n_patch,...
                (n_patch-counter) * time);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        else
           msg = sprintf(['--- Patch %i/%i --- Expected remaining ' ...
                'time: %.4f s. \r '], ...
                counter, n_patch,...
                (n_patch-counter) * time);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        
%         if i == length(x_range)
%             if j == length(y_range) 
%                 x = X(x_range(i):end,y_range(j):end,:);
%             else
%                 x = X(x_range(i):end,y_range(j):y_range(j+1)-1,:);
%             end
%         elseif j == length(y_range)
%             x = X(x_range(i):x_range(i+1)-1,y_range(j):end,:);
%         else
            x = X(x_range(i):x_range(i+1)-1,y_range(j):y_range(j+1)-1,:);
%         end
                
        % Scattering transform of the patch:
        Wop = wavelet_factory_2d(size(x), filt_opt, scat_opt);
        S_seg = scat(x, Wop);
        
        S_seg = hmm_prepare_S(S_seg, n_state);
        
        % MAP patch = ripple:
        [tmp_P_hat_ri, ~] = ...
            hmm_MAP(S_seg, theta_est_ripple, false);
        [tmp_P_hat_se, ~] = ...
            hmm_MAP(S_seg, theta_est_seabed, false);
    
        P_hat_ri = mean(mean(tmp_P_hat_ri));
        P_hat_se = mean(mean(tmp_P_hat_se));
        
        % Segmentation:
        pmap_ripple(x_range(i):x_range(i+1),y_range(j):y_range(j+1)) = P_hat_ri;
        pmap_seabed(x_range(i):x_range(i+1),y_range(j):y_range(j+1)) = P_hat_se;
        if P_hat_ri > P_hat_se
            % Ripple: %white
            segmentation(x_range(i):x_range(i+1),y_range(j):y_range(j+1)) = 1;
        else
            % Seabed: %Pink
            segmentation(x_range(i):x_range(i+1),y_range(j):y_range(j+1)) = 2;
        end
        

        
        
        counter = counter + 1;
    end
end



image_ori = figure;
image_seg = figure;
image_pRip = figure;
image_pSea = figure;
% Plot:
figure(image_ori)
imagesc(X(1:s_X(1),1:s_X(2)))
axis off
title('image')
colormap pink
drawnow

figure(image_seg)
imagesc(segmentation(1:s_X(1),1:s_X(2)))
axis off
title('segmentation')
cmap = zeros(8,3);
cmap(1:4,:) = repmat([204/255,164/255,131/255],4,1);
cmap(5:8,:) = 1;
colormap(cmap)
drawnow

%Rescalling
pmap_ripple = pmap_ripple .* 10000;
pmap_seabed = pmap_seabed .* 10000;

% Reduce max
max_seab = max(max(pmap_seabed));
pmap_seabed(pmap_seabed==max_seab) = max(max(pmap_seabed(pmap_seabed~=max_seab)));
max_rip = max(max(pmap_ripple));
pmap_ripple(pmap_ripple==max_rip) = max(max(pmap_ripple(pmap_ripple~=max_rip)));

figure(image_pRip)
imagesc(pmap_ripple(1:s_X(1),1:s_X(2)))
title('Probability of ripple')
axis off
% h = colorbar;
% set(h, 'ylim', [min(min(pmap_ripple(pmap_ripple~=0))) max(max(pmap_ripple))])
drawnow

figure(image_pSea)
imagesc(pmap_seabed(1:s_X(1),1:s_X(2)))
title('Probability of seabed')
axis off
% h = colorbar;
% set(h, 'ylim', [min(min(pmap_seabed(pmap_seabed~=0))) max(max(pmap_seabed))])
drawnow

saveas(image_ori, './Save/Plots/seg_original', 'epsc')
saveas(image_seg, './Save/Plots/seg_result', 'epsc')
saveas(image_pRip, './Save/Plots/seg_pmap_ripple', 'epsc')
saveas(image_pSea, './Save/Plots/seg_pmap_seabed', 'epsc')


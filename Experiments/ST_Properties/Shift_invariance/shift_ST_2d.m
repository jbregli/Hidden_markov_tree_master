clear all
% close all

addpath_hmt

x1 = generate_square(false, false, false);
x2 = generate_square(false, false, true);

% Precompute the wavelet transform operators that
% will be applied to the image:
% compute scattering with non-default options
filt_opt.J = 3; % scales
filt_opt.L = 2; % orientations
filt_opt.filter_type = 'morlet';
scat_opt.oversampling = 2;
scat_opt.M = 3;

[Wop, filters] = wavelet_factory_2d(size(x1), filt_opt, scat_opt);

% Call the scat function to compute the scattering
% of x using those Wop operators:
[S1, U1] = scat(x1, Wop);
[S2, U2] = scat(x2, Wop);

%% Display
%image_scat(S, true, true);


%% Data Structures
% We want to asses the persistance property over orientations 'theta'
% independently from the scales 'j'. Hence for each layer we plot a given
% orientation (same for all layers) at any available scale.

% Path to the ST
trail = {};

% find all the ST where theta = [x x x]
% All possible theta:
orientation = 1:filt_opt.L;

% Loop over the orientation:
for orient = 1:length(orientation)

    for layer=2:length(S1)
        [m,n] = size(S1{layer}.meta.theta);   
        path = orientation(orient) .* ones(m,1);
        tmp = [];

        for path_theta =1:n
            if all(S1{layer}.meta.theta(:,path_theta) == path) == 1
                tmp(end+1) = path_theta;
            end
        end
        trail{layer} = tmp;
    end

    % Plot:
    for n = 2:length(trail)
        frame_name = ['layer ' num2str(n) ' - orientation ' ...
                                            num2str(orientation(orient))];

        figure
        set(gcf,'numbertitle','off','name',frame_name) 

        for i=1:length(trail{n})
            subplot(1,length(trail{n}),i)
            imagesc(S1{n}.signal{trail{n}(i)} - S2{n}.signal{trail{n}(i)})
            title(meta2str(S1{n}.meta,trail{n}(i)));
            colormap gray; % avoids default jet colormap
        end
    end
end

figure
imagesc(x1 - x2)
colormap gray;


%%%%% TO BE INCLUDED IN FOR ALL POSSIBLE COMBINATION

% Loop on layer starting from the bottom
% pot_theta = unique(S{length(S)}.meta.theta(1,:));
% for layer=(length(S)+1) * ones(1,length(S)) - 1:length(S)

% % Display the coefficients of 2nd order with
% % j1=0, j2=2, θ1=1, θ2=5
% j1 = 0;
% j2 = 2;
% theta1 = 1;
% theta2 = 5;
% p = find( S{3}.meta.j(1,:) == j1 &...
%     S{3}.meta.j(2,:) == j2 & ...
%     S{3}.meta.theta(1,:) == theta1. & ...
%     S{3}.meta.theta(2,:) == theta2 );
% 
% imagesc(S{3}.signal{p});
% 
% % Loop on layer starting from the bottom
% for layer=(length(S)+1) * ones(1,length(S)) - [1:length(S)]
%     for j=  unique(S{layer}.meta.j(1,:))
        
    
    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script shows the 2 Population property in a Scattering Transform.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

clear all
close all

% Load an image of reasonable size:
% (e.g. 640x480)
%x = uiuc_sample;
%x = lena;

%x = generate_1st_diag(false);
x = generate_square(false, true);
%x = generate_circle(false, true);

% Precompute the wavelet transform operators that
% will be applied to the image:
% compute scattering with non-default options
filt_opt.J = 4;
filt_opt.L = 6;
filt_opt.filter_type='morlet';
scat_opt.oversampling = 2;
scat_opt.M = 2;

[Wop, filters] = wavelet_factory_2d(size(x), filt_opt, scat_opt);

% Call the scat function to compute the scattering
% of x using those Wop operators:
[Sx, Ux] = scat(x, Wop);

% % Reformat into a 3D array:
S_mat = format_scat(Sx);
 
% S_mat is a 417x60x80 matrix. The first dimension 
% in S_mat is the path index, while the second and
% third dimensions correspond to subsampled spatial
%coordinates of the original image. 

% %% Data Structures 
% 
% % Display the coefficients of 2nd order with
% % j1=0, j2=2, θ1=1, θ2=5
% j1 = 0;
% j2 = 2;
% theta1 = 1;
% theta2 = 5;
% p = find( Sx{3}.meta.j(1,:) == j1 &...
%     Sx{3}.meta.j(2,:) == j2 & ...
%     Sx{3}.meta.theta(1,:) == theta1 & ...
%     Sx{3}.meta.theta(2,:) == theta2 );
% 
% imagesc(Sx{3}.signal{p});

%% Display: 
% Display scattering coefficients
image_scat(Sx, true, true)
% opt1: scattering path display
% opt2: each node is renormalized by the norm of its parent.

% Compute the filters with 5 scales and 6 orientations
% Display all the filters
%display_filter_bank_2d(filters);

figure
imagesc(x)

%% Histogram
% for i=1:numel(Sx)
%     for j=1:numel(Sx{i}.signal)
%         figure
%         hist(reshape(Sx{i}.signal{j},numel(Sx{i}.signal{j}),1))
%     end
% end

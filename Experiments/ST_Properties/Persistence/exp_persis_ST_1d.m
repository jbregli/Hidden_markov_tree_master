%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script shows the correlation in a 1D Scattering Transform  .       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:
% Create a step-signal:
%load handel;    % loads the signal into y
y = generate_edge(false, 2^8);
N = length(y);
T = 64;       % length of the averaging window


% Precompute the wavelet transform operators that
% will be applied to the image:
% compute scattering with default options
% filt_opt = struct();
filt_opt.filter_type = 'spline_1d';
filt_opt.spline_order = 1;
% filt_opt = default_filter_options('audio', T);
% filt_opt.Q = 1;
% filt_opt.J = 1; %T_to_J(T, filt_opt);


scat_opt.M = 4;

[Wop, filters] = wavelet_factory_1d(N, filt_opt, scat_opt);

% Call the scat function to compute the scattering
% of x using those Wop operators:
[S, U] = scat(y, Wop);

% % Reformat into a 3D array:
S_mat = format_scat(S);
 

%% Display: 
S = renorm_scat(S);
S = log_scat(S);

ndisp = size(S,2);

figure

subplot(ndisp,1,1)
plot([1:size(y,1)],y)
for n = 2:ndisp
    subplot(ndisp,1,n)
    plot([1:size(S{n}.signal,2)],abs([S{n}.signal{:}]))
end
colormap gray; % avoids default jet colormap

%% Histogram
% for i=1:numel(Sx)
%     for j=1:numel(Sx{i}.signal)
%         figure
%         hist(reshape(Sx{i}.signal{j},numel(Sx{i}.signal{j}),1))
%     end
% end

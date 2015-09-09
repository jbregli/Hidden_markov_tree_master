%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         ERROR
% This script shows the shift (in)variance in a Scattering Transform.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Signal
% lenght of the signal:
len = 256;

%% ST parameters
% filt_opt = struct();
filt_opt.filter_type = 'spline_1d';
filt_opt.spline_order = 3;
% filt_opt = default_filter_options('audio', T);
% filt_opt.Q = 1;
% filt_opt.J = 1; %T_to_J(T, filt_opt);

scat_opt.M = 3;

[Wop, filters] = wavelet_factory_1d(len, filt_opt, scat_opt);

%% Shift
magnitude = {};

for i = 1:len
    % Shifted step:
    y = generate_edge(false, 2^8, i);
    N = length(y);

    % Call the scat function to compute the scattering
    % of x using those Wop operators:
    [S, U] = scat(y, Wop);

    % % Reformat into a 3D array:
    S_mat = format_scat(S);
    S = renorm_scat(S);
    S = log_scat(S);

    %Store interesting coefficients:
    for j=2:length(S)
        k = fix(i/(2^j+1)) + 1;     
        % k = fix(i/2^(length(S)-(j-1)));

        magnitude{j}(i) = S{j}.signal(k);
    end
end

%% Histogram
% for i=1:numel(Sx)
%     for j=1:numel(Sx{i}.signal)
%         figure
%         hist(reshape(Sx{i}.signal{j},numel(Sx{i}.signal{j}),1))
%     end
% end

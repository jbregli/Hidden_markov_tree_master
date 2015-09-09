%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script shows the correlation in a Dual Wavelet Tree.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:
% Create a step-signal:
%load handel;    % loads the signal into y
y = generate_edge(false, 2^8);
N = length(y);

%% DTCWT

% --- DTCWT parameters ---
biort    = 'near_sym_b_bp'; % level 1 FIR
qshift   = 'qshift_b_bp';   % Q-shift FIR


% meanY   = mean(y(:));
% y       = y-meany;
N_levels = 5;
[Yl, Yh, Yscale] = wavexfm(y, N_levels, biort);


% Display:
ndisp = size(Yh,1);

figure
set(gcf,'numbertitle','off','name','DWT') 

subplot(ndisp+1,1,1)
plot([1:size(y,1)],y)
for n = 0:ndisp-1
        subplot(ndisp+1,1,n+2)
        plot([1:size(Yh{ndisp-n},1)],[abs(Yh{ndisp-n})])
end


%% SOLUTION 2:
% J = 3;
% dwt1 = dddtree('dwt',y,J,'db8');
% 
% ndisp = size(dwt1.cfs,2);
% 
% figure
% set(gcf,'numbertitle','off','name','DWT') 
% 
% subplot(ndisp+1,1,1)
% plot([1:size(y,1)],y)
% for n = 0:ndisp-1
%         subplot(ndisp+1,1,n+2)
%         plot([1:size(dwt1.cfs{ndisp-n},1)],[abs(dwt1.cfs{ndisp-n})])
% end

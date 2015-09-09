%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script shows the correlation in a Dual Wavelet Complex Tree.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Initialization:
% Create a step-signal:
y = generate_edge(false, 2^8);
N = length(y);

%% DTCWT

% --- DTCWT parameters ---
biort    = 'near_sym_b_bp'; % level 1 FIR
qshift   = 'qshift_b_bp';   % Q-shift FIR


% meanY   = mean(y(:));
% y       = y-meany;
N_levels = 6;
[Yl, Yh, Yscale] = dtwavexfm(y, N_levels, biort, qshift);


%Display:
ndisp = size(Yh,1);

figure;
set(gcf,'numbertitle','off','name','DWCT') 

subplot(ndisp+1,1,1);
plot([1:size(y,1)],y);
for n = 0:ndisp-1
        subplot(ndisp+1,1,n+2);
        plot([1:size(Yh{ndisp-n},1)],[abs(Yh{ndisp-n})]);
end
title('DWCT')


%% SOLUTION 2:
% J = 3;
% dwt1 = dddtree('cplxdt',y,J,'dtf3');
% 
% 
% ndisp = size(dwt1.cfs,2);
% 
% figure
% set(gcf,'numbertitle','off','name','DWCT') 
% 
% subplot(ndisp+1,1,1)
% plot([1:size(y,1)],y)
% for n = 0:ndisp-1
%         signal = dwt1.cfs{ndisp-n}(:,:,1) + 1j *  dwt1.cfs{ndisp-n}(:,:,2);
%         
%         subplot(ndisp+1,1,n+2)
%         plot([1:size(dwt1.cfs{ndisp-n},1)],[abs(signal)])
% end
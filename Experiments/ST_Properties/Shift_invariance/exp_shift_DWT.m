%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         OK
% This script shows the shift (in)variance in a Dual Wavelet Tree.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Signal
% lenght of the signal:
len = 180;

%% DTCWT parameters ---
biort    = 'near_sym_b_bp'; % level 1 FIR
qshift   = 'qshift_b_bp';   % Q-shift FIR

N_levels = 4;

%% Shift
magnitude = {};

for i = 1:len
    
    % Shifted step:
    y = generate_edge(false, 2^8, i);
    N = length(y);

    % DTW transform
    [Yl, Yh, Yscale] = wavexfm(y, N_levels, biort);

    %Store interesting coefficients:
    for j=1:length(Yh)
        k = fix(i/(2^j+1)) + 1;  
        magnitude{j}(i) = Yh{j}(k);
    end
    
    % Ignore full 1 and full 0
    if i < 6 ||i > len -6
       for j=1:length(magnitude)
           magnitude{j}(i) = 0;
       end
    end
        
end

%% Display:
ndisp = length(magnitude);

figure;
set(gcf,'numbertitle','off','name','DWT - Shift') 
for n = 1:ndisp
        subplot(ndisp,1,n);
        plot([1:size(magnitude{n},2)],[magnitude{ndisp-(n-1)}]);
end
title('DWCT')


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

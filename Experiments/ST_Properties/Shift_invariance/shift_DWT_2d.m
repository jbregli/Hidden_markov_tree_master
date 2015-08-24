clear all
close all

addpath_hmt

x1 = generate_square(false, false, false);
x2 = generate_square(false, false, true);

% --- DTCWT parameters ---
biort    = 'near_sym_b_bp'; % level 1 FIR
qshift   = 'qshift_b_bp';   % Q-shift FIR


% meanY   = mean(y(:));
% y       = y-meany;
N_levels = 5;
[Yl1, Yh1, Yscale1] = wavexfm(x1, N_levels, biort);
[Yl2, Yh2, Yscale2] = wavexfm(x2, N_levels, biort);

%% Display:
ndisp = size(Yh1,1);

figure
set(gcf,'numbertitle','off','name','DWT') 

subplot(ndisp+1,1,1)
imagesc(x1 - x2)
for n = 0:ndisp-1
        subplot(ndisp+1,1,n+2)
        imagesc(abs(Yh1{ndisp-n}) - abs(Yh2{ndisp-n}))
end
colormap gray
        
    
    




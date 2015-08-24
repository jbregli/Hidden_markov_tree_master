path_to_hmt = fileparts(mfilename('fullpath'));

addpath(genpath('./'));

addpath(fullfile(path_to_hmt, '../../scatnet-master'));
addpath_scatnet;

clear path_to_hmt;
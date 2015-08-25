path_to_hmt = fileparts(mfilename('fullpath'));

addpath(genpath('./'));

path_to_scatnet = fullfile(path_to_hmt, '../../scatnet-master');

if exist(path_to_scatnet,'dir')
    addpath(fullfile(path_to_hmt, '../../scatnet-master'));
    addpath_scatnet;
end

clear path_to_hmt;
clear path_to_scatnet;
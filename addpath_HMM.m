path_to_hmt = fileparts(mfilename('fullpath'));

addpath(fullfile(path_to_hmt, 'Generate'));

addpath(fullfile(path_to_hmt, 'Properties'));
% addpath(fullfile(path_to_hmt, 'Properties/Persistence'));
addpath(fullfile(path_to_hmt, 'Properties/Shift'));
addpath(fullfile(path_to_hmt, 'Properties/TwoPop'));
addpath(fullfile(path_to_hmt, 'Properties/Distribution'));
addpath(fullfile(path_to_hmt, 'Properties/Correlation'));
addpath(fullfile(path_to_hmt, 'Properties/Histogram'));

addpath(fullfile(path_to_hmt, 'BayesNet'));

addpath(fullfile(path_to_hmt, 'Hmm'));
addpath(fullfile(path_to_hmt, 'Hmm/E_step'));
addpath(fullfile(path_to_hmt, 'Hmm/E_step/Beta'));
addpath(fullfile(path_to_hmt, 'Hmm/E_step/CP'));
addpath(fullfile(path_to_hmt, 'Hmm/M_step'));
addpath(fullfile(path_to_hmt, 'Hmm/Prepare'));

addpath(fullfile(path_to_hmt, 'Hmm2'));

addpath(fullfile(path_to_hmt, 'Script'));
addpath(fullfile(path_to_hmt, 'Test'));

addpath(fullfile(path_to_hmt, '../scatnet-master'));
addpath_scatnet;


clear path_to_hmt;


function [data, group ] = load_svm(path_to_set, class)
% load_svm: LOAD DATA OF A CLASS TO MEET 'SVMTRAIN' REQUIREMENTS
%
%   --------
%   INPUTS:
%   --------
%   - path_to_set: string or cell{string, int}
%       Two types of input are accepted. Either a string giving a path to a
%       folder containing a set of similar images. Or a cell giving the
%       name of a known generator ('square' or 'circle') and the number of
%       images to generate.
%   - class: int
%       Reference of the class
%
%   --------
%   OUTPUTS:
%   --------
%   - data: matrix
%       Matrix of training data, where each row corresponds to an 
%       observation or replicate, and each column corresponds to a feature
%       or variable. 
%   - group: vector
%       Each element of Group specifies the group of the corresponding row
%       of Training. Group should divide Training into two groups. 
%       Group has the same number of elements as there are rows in 
%       Training. 
%

    %% Initialization:
    % Display
    fprintf(' * Loading data: \n')

    % List all the images:
    allFiles = dir(fullfile(path_to_set, '*.png'));
    allNames = {allFiles.name};

    % Input caracteristics:
    n_image = length(allNames);
    s_image = size(imread(fullfile(path_to_set, allNames{1})));
    
    % Data & group matrix:
    data = zeros(n_image, prod(s_image));
    group = zeros(n_image, 1);

    for i = 1:n_image
        data(i,:) = ...
            reshape(...
                im2double(imread(fullfile(path_to_set, allNames{i}))),...
                prod(s_image),1);

        group(i,1) = class;
    end  
end


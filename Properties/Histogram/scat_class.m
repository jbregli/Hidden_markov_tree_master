function [ transform ] = scat_class(path_to_set, filt_opt, scat_opt)

% scat_class: COMPUTES THE SCATTERING TRANSFORM OF A BATCH OF SIMILAR
%             IMAGES
%
%   Given a set of images (path to it or generator) computes the STs.
%
%   --------
%   INPUTS:
%   --------
%   - path_to_set: string or cell{string, int}
%       Two types of input are accepted. Either a string giving a path to a
%       folder containing a set of similar images. Or a cell giving the
%       name of a known generator ('square' or 'circle') and the number of
%       images to generate.
%   - filt_opt: (optional) struct
%       Filtering options for the scattering transform
%   - scat_opt: (optional) struct)
%       Scattering options for the scattering transform
%
%   --------
%   OUTPUTS:
%   --------
%   - transform: cell
%       Cell of the same length as the number of inputs storing the
%       scattering transform
%
%   WARNING: RAM usage (?) to avoid potential issues work with small
%            (100/200) batches.

    %% Initialization:
    if nargin < 2
        filt_opt = struct();
    end
    if nargin< 3
        scat_opt = struct();
    end
    
    % Display
    reverseStr = '';
    fprintf(' * Image processing: \n')

    %% Case where the input is a path:
    if ischar(path_to_set)

        % List all the images:
        allFiles = dir(fullfile(path_to_set, '*.png'));
        allNames = {allFiles.name};

        % Scattering transform:
        % Initilization w/ first image:
        x = im2double(imread(fullfile(path_to_set, allNames{1})));
        % Pre-compute the WT op that will be applied to the image:
        Wop = wavelet_factory_2d(size(x), filt_opt, scat_opt);

        tic;
        S = scat(x, Wop);
        time = toc;

        % Stock STs in a cell:
        transform = [{} {S}];
        
        % LOOP OVER THE IMAGES:
        for i = 2:length(allNames)
            % Print time remaining:
            msg = sprintf('--- Image %i/%i --- Expected remaining time: %.4f s. \r ' ,...
                i, length(allNames), (length(allNames)-i) * time);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            
            % ST:
            [S, U] = scat(im2double(imread(fullfile(path_to_set, allNames{i}))), Wop);
            
            transform = [transform {S}];
        end

    %% Case where the input is a cell:
    elseif iscell(path_to_set)
        if length(path_to_set) < 3
            size_im = [640 640];
        else
            size_im = path_to_set{3};
        end
        
        if strcmp(path_to_set{1},'square')
            % Generate path_to_set{2} number of square and their ST
            Wop = wavelet_factory_2d(size_im, filt_opt, scat_opt);

            x = generate_square(true, true, true, false, size_im);
            % empty, noise, translate, rotate, size)

            tic;
            S = scat(x, Wop);
            time = toc;

            % Stock STs in a cell:
            transform = [{} {S}];
            % GENERATE THE REQUIERED NUMBER OF IMAGES:
            for i = 2:path_to_set{2}
                msg = sprintf('--- Image %i/%i --- Expected remaining time: %.4f s. \r ' ,...
                    i, path_to_set{2}, (path_to_set{2}-(i-1)) * time);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                
                % ST:
                x = generate_square(true, true, true, false, size_im);
                [S, U] = scat(x, Wop);

                transform = [transform {S}];
            end
        else
            disp('Generator not implemented yet')
        end

    %% Error catching:
    else
        disp('Invalid input')
    end
    
    
end


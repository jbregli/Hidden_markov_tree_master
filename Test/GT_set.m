function [set_S] = GT_set(S, theta_gt, n_image)
% GT_theta:  CREATE A SET OF GENERATE IMAGES
%   Given a scattering transform (cell of structure), this functions 
%   associates to each node the required fields for the HMM modeling.
%
%   --------
%   INPUTS:
%   --------
%   - S: cell(struct)
%       Set of structures obtained with the function 'scat' of the 
%       'scatnet' lib.     
%   - theta_gt: cell(struct)
%       Structure with the same organisation as 'S' the 'scatnet' lib
%       containing the model parameters.
%   - n_image: int
%       Number of image to be generated.
%   
%   --------
%   OUTPUTS:
%   --------
%  - set_S: cell{cell(struct)}
%       Set of structures containing sampled images.    #
%
%   --------
%   IMPROVEMENTS:
%   --------
%

    %% Initialization:   
           
    %% Generate the images:
    set_S = cell(1, n_image);
    for i=1:n_image
        set_S{i} = GT_tree(S, theta_gt);
    end
            
            
    
    

end


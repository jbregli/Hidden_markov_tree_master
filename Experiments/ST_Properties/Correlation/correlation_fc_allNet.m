function [corr_strct] = correlation_fc_allNet(S, sct_plot, verbose)
% correlation_fc_allNet: CALL 'correlation_fc' FOR ALL THE NODE OF THE
%                        NETWORK
%   Given a scattering transform, this function computes the correlation
%   coefficient of all the father nodes with their respective children.
%
%   --------
%   INPUTS:
%   --------
%   - S: struct
%       Structure obtained with the function 'scat' of the 'scatnet' lib.
%
%   --------
%   OUTPUTS:
%   --------
%   - corr_strct: struct
%       Structure holding the correlation score between the father node's
%       scattering coefficient and those of each of its child. if given as
%       an input the old structure is updated, otherwise it's created.

    %% Initialization:
    if nargin < 2
        sct_plot = false;
    end
    if nargin < 3
        verbose = false;
    end
   
    corr_strct = {};
    
    
    %% Loop over the network
    for f_layer=3:(length(S)-1)
        % Loop over the nodes of a layer
        for f_index=1:length(S{f_layer}.signal)
            % Call 'correlation_fc (no plots)         
            corr_strct = correlation_fc(S, f_layer, f_index, sct_plot, ...
                                                      corr_strct, verbose);
                                                        
        end
    end  
end


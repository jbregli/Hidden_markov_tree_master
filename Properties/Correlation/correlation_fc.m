function [corr_strct] = correlation_fc(S, f_layer, f_index, ...
                                        sct_plot, corr_strct, verbose)
% correlation_fc: CORRELATION BETWEEN THE SCATTERING COEFFICIENTS OF ALL
%                 THE CHILDREN OF A NODE
%
%   Given a layer and an index for the father node this function compute
%   the correlation between its scattering coefficient and those of each of
%   its child. 
%
%   --------
%   INPUTS:
%   --------
%   - S: struct
%       Structure obtained with the function 'scat' of the 'scatnet' lib.
%   - f_layer: int
%       Layer of the father node
%   - f_index: int
%       Index of the father node within the layer 'f_layer'
%   - sct_plot: (optional) bool
%       Do you want any figure to be plotted?
%   - corr_strct: (optional) struct
%       Structure holding the correlation score between the father node's
%       scattering coefficient and those of each of its child.
%
%   --------
%   OUTPUTS:
%   --------
%   - corr_strct: struct
%       Structure holding the correlation score between the father node's
%       scattering coefficient and those of each of its child. if given as
%       an input the old structure is updated, otherwise it's created.

    %% Initialization:
    % Structure saving the correlation scores
    if nargin < 4
        sct_plot = true;
    end
    
    if nargin < 5
        corr_strct = {};
    end
    if not(isfield(corr_strct, 'layer'))
        corr_strct.layer = {}; 
    end
    if not(isfield(corr_strct, 'scale'))
        corr_strct.scale = {}; 
    end
    if not(isfield(corr_strct, 'theta'))
        corr_strct.theta = {}; 
    end
    if not(isfield(corr_strct, 'corr'))
        corr_strct.corr = {}; 
    end
    
    if nargin < 6
        verbose = true;
    end
    
    %% Find path to the father:
    f_scale = S{f_layer}.meta.j(:,f_index)';
    f_theta = S{f_layer}.meta.theta(:,f_index);

    %% Find children:
    f_l_scale = length(S{f_layer}.meta.j(:,f_index));
    f_l_theta = length(S{f_layer}.meta.theta(:,f_index));
    c_l_scale = length(S{f_layer+1}.meta.j(:,1));    
    c_l_theta = length(S{f_layer+1}.meta.theta(:,1));

    % Find in the next layer the index of the children of the father node:
    % Same scale path:
    if f_layer == 1
        % all the nodes of the 2nd layers are children:
        c_index = 1:length(S{f_layer+1}.meta.j);
    else
        temp = ismember(S{f_layer+1}.meta.j((c_l_scale-(f_l_scale-1)):c_l_scale,:)',...
                                                          f_scale, 'rows');
        index_scale = find(temp == 1);
        % Same orientation
        temp = ismember(S{f_layer+1}.meta.theta((c_l_theta-(f_l_theta-1)):c_l_theta,:)',...
                                                          f_theta', 'rows');
        index_theta= find(temp == 1);
        % Find the matching index:
        c_index = intersect(index_scale, index_theta);
    end

    %% Case where this node has no children:
    if isempty(c_index)
        if verbose == true
            disp(['The node ' num2str(f_index) ' at layer ' num2str(f_layer) ...
                  ' does not have any child.'])
        end
    else
        %% Correlation score and scatter plot:
        if sct_plot == true
            % Scattering transform plot init:
            figure(1)
            
            %Scatter plot init:
            figure(2)
            xlabel('ST coefficient - Father');
            ylabel('ST coefficient - Child');
        end        
        
        % X = father:
        X = reshape(S{f_layer}.signal{f_index},...
                                      numel(S{f_layer}.signal{f_index}),1);       
        for i=1:length(c_index)
            % Y = child:
            Y = reshape(S{f_layer+1}.signal{c_index(i)},...
                                 numel(S{f_layer+1}.signal{c_index(i)}),1);

            % Correlation:
            correlation = corr(X, Y);
                   
            % Plot:
            if sct_plot == true
                if i > 1
                    pause(2);
                end
                figure(1)
                % Father:
                subplot(3,1,1)
                imagesc(S{f_layer}.signal{f_index})
                title('Father');
                colormap gray;
                 % Child:
                subplot(2,1,2)
                imagesc(S{f_layer+1}.signal{c_index(i)})
                title('Son');  
                colormap gray;
                
                figure(2)
                frame_name = ['Correlation layers ' mat2str([f_layer f_layer+1])  ...
                              ' - orientation ' mat2str(S{f_layer+1}.meta.theta(:,c_index(i)))...
                              ' - scale ' mat2str(S{f_layer+1}.meta.j(:,c_index(i)))];

                set(gcf,'numbertitle','off','name',frame_name)     
                scatter(X,Y)
                xlabel('Father ST coef')
                ylabel('Child ST coef')
                title(['Correlation = ' num2str(correlation)]);
                
                
                
                
            end
            
            % Update the structure:
            corr_strct.layer= [corr_strct.layer {f_layer+1}];
            corr_strct.scale = [corr_strct.scale {S{f_layer+1}.meta.j(:,c_index(i))}];
            corr_strct.theta = [corr_strct.theta {S{f_layer+1}.meta.theta(:,c_index(i))}];
            corr_strct.corr = [corr_strct.corr {correlation}];
        end
    end
end


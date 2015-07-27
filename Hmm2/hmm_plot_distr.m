function hmm_plot_distr(set_S, theta_est, cv_achieved, layer, scale, x, y, fig)
% hmm_plot_distr: PLOT THE ESTIMATED PARAMETER AGAINST THE GROUND TRUTH
%   If the ground truth is given. Otherwise it just plots the estimated 
%   parameters.
%
%   --------
%   INPUTS:
%   --------
%   - theta_est: cell(struct)
%       Cell of structures with the same organisation as 'S' the 'scatnet' 
%       lib where the estimated paramters are stored.
%           .proba, .epsilon, . mu, .sigma, .distr
%   - layer: (optional) int (default: random)
%       Layer number of 'var_TBC'. Help for debugging. Display only if 
%       'verbose' is true.
%   - scale: (optional) int (default: random)
%       Scale number of 'var_TBC'. Help for debugging. Display only if 
%       'verbose' is true.
%
%   --------
%   OUTPUTS:
%   --------
%
%   --------
%   ISSUES:
%   --------
%
%   --------
%   TODO:
%   --------
%
    %% Preparation:   
    close all
    
    % Arguments (part 1):
    if ~exist('fig','var')
        fig = figure;
    end

    % Sizes:
    n_image = length(set_S);
    n_layer = length(theta_est);
    n_state = size(theta_est{1}.proba{1}, 3);
    s_image = size(theta_est{1}.proba{1}(:,:,1));
    
    n_scale = zeros(1,n_layer);
 
    for lay=1:n_layer
        n_scale(1,lay) = length(theta_est{lay}.proba);
    end
    
    % Arguments (part 2):
    if ~exist('layer','var')
        layer = randi(n_layer);
    end
    if ~exist('scale','var')
        scale = randi(n_scale(layer));
    end
    
    if ~exist('x','var')
        % Select a random pixel in the image:    
        x = randi(s_image(1));
    end  
    if ~exist('y','var')
        % Select a random pixel in the image:    
        y = randi(s_image(2));
    end
    
    % List of colors available:
    colors = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];
    
    %% Plot:
    % Title:
    p_title = sprintf('Distribution at layer %i and scale %i', layer, scale);
               
    % Means and variances for the estimate:
    if layer == 1
        mu_est = squeeze(theta_est{layer}.proba{scale}(x,y,:));
        sigma_est = [0 0];
    else
        mu_est = squeeze(theta_est{layer}.mu{scale}(x,y,:));
        sigma_est = squeeze(theta_est{layer}.sigma{scale}(x,y,:));
    end
        
    for im = 1:n_image
        cluster(im) = set_S{im}{layer}.signal{scale}(x,y);
    end
       
    figure(fig)
    subplot(2,1,1);
    %hax = get(gca,'axes'); 
       
    hold on
    hist(cluster)
    
    for state=1:n_state
        % Mean bar:
        line([mu_est(state) mu_est(state)], ...
            0.75*get(gca,'YLim'), 'Color', colors(state),...
            'LineWidth',5);
        % Variance bar:
        line([mu_est(state)-sigma_est(state) mu_est(state)+sigma_est(state)], ...
            [1 1], ...
            'Color', colors(state), ...
            'LineWidth',5);
    end

    title(sprintf('Pixel (%i,%i) - Layer %i - Scale %i', x,y, layer, scale));
    hold off
    
    subplot(2,1,2);
    imagesc(squeeze(cv_achieved{layer}.epsilon{scale}(x,y,:,:)))
    colormap gray;
    
    drawnow
       
end


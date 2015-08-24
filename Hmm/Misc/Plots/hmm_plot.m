function hmm_plot(theta_est, theta_gt, parameter, layer, scale)
% hmm_plot_vsGT: PLOT THE ESTIMATED PARAMETER AGAINST THE GROUND TRUTH
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
%   - theta_gt: (optional) cell(struct) (default: false)
%       Same architecture as 'theta_est' but storing the ground truth 
%       values of the parameters.
%   - parameter: (optional) string (default: 'distribution')
%       Function name from where the checking is done. Help for debugging.
%       Display only if 'verbose' is true.
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
    close all;
    
    % Arguments (part 1):
    if ~exist('theta_gt','var')
        theta_gt = {};
    end
    if ~exist('parameter','var')
        parameter = 'distribution';
    end
    
    % Sizes:
    n_layer = length(theta_est);
    n_state = size(theta_est{1}.proba{1}, 3);
    s_image = size(theta_est{1}.proba{1}(:,:,1));
    
    n_scale = zeros(1,n_layer);
    

    for layer=1:n_layer
        n_scale(1,layer) = length(theta_est{layer}.proba);
    end
    
    % Arguments (part 2):
    if ~exist('layer','var')
        layer = randi(n_layer);
    end
    if ~exist('scale','var')
        scale = randi(n_scale(layer));
    end

    %% Plot:
    % Title:
    p_title = sprintf('%s at layer %i and scale %i', parameter, layer, scale);
    
    % Plot:
    switch parameter
        % Distribution:
        case 'distribution'
            distrib = theta_est{layer}.distr;
            
            % Select a random pixel in the image:
            x = randi(s_image(1));
            y = randi(s_image(2));
            
            if ~isempty(theta_gt)
                % Means and variances for the estimate:
                mu_est = squeeze(theta_est{layer}.mu{scale}(x,y,:));
                sigma_est = squeeze(theta_est{layer}.sigma{scale}(x,y,:));
            else
                % Means and variances for the estimate:
                mu_est = squeeze(theta_est{layer}.mu{scale}(x,y,:));
                sigma_est = squeeze(theta_est{layer}.sigma{scale}(x,y,:));
                % Means and variances for the ground truth:
                mu_gt = squeeze(theta_gt{layer}.mu{scale}(x,y,:));
                sigma_gt = squeeze(theta_gt{layer}.sigma{scale}(x,y,:));
            end

            switch distrib
                case 'MixtGauss'
                    if ~isempty(theta_gt)
                        Xmin = mu_est - 3 * sigma_est;
                        Xmax = mu_est + 3 * sigma_est;

                        figure('name',p_title);
                        for m=1:n_state
                            X = Xmin(m):0.1:Xmax(m);

                            estimate = normpdf(X, mu_est(m), sigma_est(m));

                            subplot(1,n_state,m)
                            plot(X, estimate, 'blue')
                            title(sprintf('State %i at pixel (%i,%i)',m,x,y));
                        end
                    else
                        Xmin = min(mu_est - 3 * sigma_est, mu_gt - 3 * sigma_gt);
                        Xmax = max(mu_est + 3 * sigma_est, mu_gt + 3 * sigma_gt);

                        figure('name',p_title);
                        for m=1:n_state
                            X = Xmin(m):0.1:Xmax(m);

                            estimate = normpdf(X, mu_est(m), sigma_est(m));
                            ground_T = normpdf(X, mu_gt(m), sigma_gt(m));

                            subplot(1,n_state,m)
                            plot(X, estimate, 'blue', ...
                                 X, ground_T, 'red')
                            title(sprintf('State %i at pixel (%i,%i)',m,x,y));
                        end    
                    end
                    
                case 'FoldedGauss'
                    disp('TBD')
            end
        case 'TBD'
            disp('TBD')           
    end
end


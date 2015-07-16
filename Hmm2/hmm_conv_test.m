function cv_achieved = hmm_conv_test(theta, theta_old, step, ...
                                         mixing, sensibility)
% hmm_conv_test: CONVERGENCE TEST FOR EM ALGORITHM
%   Given a theta at step n and theta at step n-1, this function asseses if
%   convergence has occcured.

    %% Preparation:
    % optional variables:
    if ~exist('mixing','var')
        mixing = 10;
    end
    if ~exist('sensibility','var')
        sensibility = 1e-4;
    end

    % Sizes:
    n_layer = length(theta);

    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(theta{layer}.proba);
    end

    % Message:
    msg = '';

    %% Convergence test:
    % If the number of step is too small then cv didn't occur (burning
    % time)
    if step < mixing
        cv_achieved = false;
    else
        layer = 1;
        cv_achieved = true;
        % Loop over the layers:
        while (layer <= n_layer && cv_achieved)
            %  Loop over the scales at 'layer':
            scale = 1;
            while (scale <= n_elmt(1,layer) && cv_achieved)
                % Delta proba:
                delta.proba = theta{layer}.proba{scale} - theta_old{layer}.proba{scale};
                % Delta epsilon:
                delta.epsilon = theta{layer}.epsilon{scale} - theta_old{layer}.epsilon{scale};
                % Delta mu:
                delta.mu = theta{layer}.mu{scale} - theta_old{layer}.mu{scale};
                % Delta sigma:
                delta.sigma = theta{layer}.sigma{scale} - theta_old{layer}.sigma{scale};

                i = 1;
                fields = fieldnames(delta);

                while (i < numel(fields) && cv_achieved == true)
                    if any( ...
                            delta.(fields{i})(:) > ...
                            sensibility * ones(size(delta.(fields{i})(:))))
                        cv_achieved = false;
                    end
                    i = i+1;
                end
                scale = scale +1;
            end
            layer = layer+1;
        end
    end
end


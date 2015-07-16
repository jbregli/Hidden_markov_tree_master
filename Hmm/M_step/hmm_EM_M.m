function [set_S, theta] = hmm_EM_M(set_S, theta)
% hmm_EM_M: COMPUTE THE MAXIMISATION STEP OF THE EXPECTATION MAXIMISATION
%           ALGORITHM
%   
%   See "Wavelet-based statistical signal processing using hidden markov
%   models" for more details on the algorithm
%
%   --------
%   INPUTS:
%   --------
%   - set_S: cell{cell(struct)}
%       Set of structures obtained with the function 'hmm_EM_E'.      
%
%   --------
%   OUTPUTS:
%   --------
%   - set_S: cell(struct)
%       Set of structures updated
%
%   --------
%   IMPROVEMENTS:
%   --------
%   - Lot of optimization possible in the computation of the betas
%   (beta_exclude node is computed several times, useless for loops (?),
%   useless variables...

     %% Initialization:
     % Sizes:
    n_image = length(set_S);
    n_layer = length(set_S{1});
    n_elmt = zeros(1,n_layer);
    for layer=1:n_layer
        n_elmt(1,layer) = length(set_S{1}{layer}.signal);
    end
    n_state = size(theta{1}.proba{1}, 3);
    
    % +++ Sanity check:
    Stheta = {'proba' 'epsilon' 'mu' 'sigma'}; 
    for layer=1:length(theta)
        for i=1:numel(Stheta)
            for scale=1:length(theta{layer}.(Stheta{i}))
                if max(max(max(isnan(theta{layer}.(Stheta{i}){scale}))))
                    disp(['EM_M: init check NAN theta_' Stheta{i} ' at layer ' num2str(layer) ...
                          ' and scale ' num2str(scale)])
                end
                if max(max(max(theta{layer}.(Stheta{i}){scale} == 0))) ==1
                    disp(['EM_M: init check 0 theta_' Stheta{i} ' at layer ' num2str(layer) ...
                          ' and scale ' num2str(scale)])
                end                
                
            end
        end
    end
    % +++ ---> issue with exclude node
    

    %% Update the parameter vector theta:
    for layer=1:n_layer        
        for scale=1:n_elmt(1,layer)
            % Summing temporary variables:
            tmp_proba = zeros(size(theta{layer}.proba{1}));
            tmp_epsilon = zeros(size(theta{layer}.epsilon{1}));
            tmp_mu = zeros(size(theta{layer}.mu{1}));
            tmp_sigma = zeros(size(theta{layer}.sigma{1}));
            
            % Resizing temporary variables:
            tmp_father = zeros(size(theta{layer}.epsilon{1}));
            tmp_W = zeros(size(theta{layer}.mu{1}));
            
            % PROBA OF STATE:
            % p_{S_i}(m) :
            for im=1:n_image
                tmp_proba = tmp_proba + ...
                              set_S{im}{layer}.hmm{scale}.condProb.givenVis;
            end
            
            % Normalize 'proba':
            tmp_proba = tmp_proba / n_image;
            
            % Try to avoid overflow by tresholding:
            tmp_proba = tmp_proba.*(tmp_proba>1e-4) + 1e-4*(tmp_proba<=1e-4);

            % Proba:
            theta{layer}.proba{scale} = tmp_proba; 

            % +++ Sanity check:
            if max(max(max(isnan(tmp_proba))))
                disp(['EM_M: NAN tmp_proba at layer ' num2str(layer) ...
                      ' and scale ' num2str(scale)])
            end
            if max(max(max(tmp_proba == 0))) ==1
                disp(['EM_M: 0 tmp_proba at layer ' num2str(layer) ...
                      ' and scale ' num2str(scale)])
            end
            if max(max(max(abs(tmp_proba) == inf))) ==1
                disp(['EM_M: inf tmp_proba 1 at layer ' num2str(layer) ...
                      ' and scale ' num2str(scale)])
            end
             % +++
                                
            % The means and variances of the mixed gaussians as well as the
            % probabilities of transition are not define for the root
            % layer:
            if layer > 1                               
                % PROBA OF TRANSITION & MEAN & VARIANCE:
                % epsilon_{i,rho(i)}^{mn} & mu_{i} & sigma_{i} :
                for im=1:n_image
                    % Sum epsilon:
                    tmp_epsilon = tmp_epsilon + ...
                        set_S{im}{layer}.hmm{scale}.condProb.givenVisAndParents;                  
                    
                    % Resizing the ST coefficient for matrix product:                   
                    tmp_W = repmat(set_S{im}{layer}.signal{scale},1,1,n_state);

                    % Sum mean:
                    tmp_mu = tmp_mu + tmp_W ...
                        .* set_S{im}{layer}.hmm{scale}.condProb.givenVis;
                    
                    % Sum variance:
                    tmp_sigma = tmp_sigma ...
                       + ((tmp_W - theta{layer}.mu{scale}).^2 ...
                        .* set_S{im}{layer}.hmm{scale}.condProb.givenVis);  
                end
                
                % Resizing the father probability for matrix product:
                % from [s_im 1] to [s_im n_state, n_state]
                for m=1:n_state
                    for n=1:n_state  
                        tmp_father(:,:,n,m) = theta{layer-1}.proba{set_S{1}{layer}.hmm{scale}.parent}(:,:,n);
                    end
                end
                
                % Normalize 'epsilon':
                tmp_epsilon = tmp_epsilon ./ (n_image .* tmp_father);               
                % +++ Ensure everything sum to 1:
                tmp_epsilon = hmm_prepare_normalEps(tmp_epsilon);
                
                % Normalize 'mu': 
                tmp_mu = tmp_mu ./ (n_image * theta{layer}.proba{scale});
                % Try to avoid overflow by tresholding:
                %tmp_mu = tmp_mu.*(tmp_mu>1e-4) + 1e-4*(tmp_mu<=1e-4);
                
                % Normalise 'sigma':
                tmp_sigma = tmp_sigma ./ (n_image * theta{layer}.proba{scale});
                % Try to avoid overflow by tresholding:
                %tmp_sigma = tmp_sigma.*(tmp_sigma>1e-4) + 1e-4*(tmp_sigma<=1e-4);                          
                
                %% Updates:   
                % Epsilon: 
                theta{layer}.epsilon{scale} = tmp_epsilon;
                % Mu:
                theta{layer}.mu{scale} = tmp_mu; 
                % Sigma:
                theta{layer}.sigma{scale} = tmp_sigma;
            end
            % +++ Try with non sequential updates ---> moved before the
            % 'if'
%             % Proba:
%             theta{layer}.proba{scale} = tmp_proba;
            
            %% +++ Sanity check:
            Stheta = fieldnames(theta{layer});
            for i=1:numel(Stheta)
                if not(strcmp('distr',Stheta{i}))
                    if max(max(max(isnan(theta{layer}.(Stheta{i}){scale}))))
                        disp(['EM_M: NAN theta_' Stheta{i} ' at layer ' num2str(layer) ...
                            ' and scale ' num2str(scale) ...
                            ' for image ' num2str(im)])
                    end
                    if max(max(max(theta{layer}.(Stheta{i}){scale} ==0))) == 1
                        disp(['EM_M: 0 theta_' Stheta{i} ' at layer ' num2str(layer) ...
                            ' and scale ' num2str(scale) ...
                            ' for image ' num2str(im)])
                    end
                end
            end
            % +++ ---> issue with sigma ans mu
        end
    end
            
    

end


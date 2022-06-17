function [pca_out, pca_array] = snr_pca(snr_by_bhv, nv)
% Accepts as input the snr_out_by_bhv output of snr_all_subj_sess.m and
% performs PCA on a behaviors x regions array to search for
% cananonical behavioral circuits, canoncical neural reigonal circuits, and
% the relationship between the two sets.
%
% PARAMETERS
% ----------
% snr_by_bhv : snr_out_by_bhv output of snr_all_subj_sess.m
% Name-Value Pairs (nv)
%   'obsevations' : ('single_session' | 'session_mean'), specify what
%                   would constitute an observation in a PCA model. If
%                   'single_session', SNR values from each behavior in each
%                   session are used. If 'session_mean', SNR values are
%                   first averaged across sessions within each behavior.
%   'normalize'   : (1 | 0 (default)), specify whether or not to z-score
%                   columns of input array to PCA. If true, this is
%                   equivalent to correlation matrix PCA.
%
% RETURNS
% -------
% pca_out   : PCA results; output of pca_svd.m function.
% pca_array : Array passed to PCA.
% 
% Author: Jonathan Chien 6/2/22.


arguments
    snr_by_bhv
    nv.observations string = 'single_session' 
    nv.normalize = false
end

% Get names of behaviors and array sizes. Initialize array that will be
% input to PCA/SVD.
behaviors = fieldnames(snr_by_bhv);
n_behaviors = length(behaviors);
pca_array = [];

% For each behavior/session, get observed snr metric, as well as mean and
% std of asociated null distribution.
for i_bhv = 1:n_behaviors
    % Get field of input struct corresponding to current behavior.
    curr_bhv = horzcat(snr_by_bhv.(behaviors{i_bhv}).by_region)';
    
    % Optional averaging step.
    if strcmp(nv.observations, 'session_mean')
        pca_array = [pca_array; mean(curr_bhv, 1, 'omitnan')];

    % Eschew averaging and use individual sessions as input to PCA.
    elseif strcmp(nv.observations, 'single_session')
        % Check for NaNs. If found, issue warning and remove them.
        nan_obs = any(isnan(curr_bhv), 2);
        if any(nan_obs) 
            curr_bhv(nan_obs,:) = [];  
            warning("Z scores for all regions from session(s) %s from " + ...
                    "%s removed due to presence of NaN values.", ...
                    strjoin(string(find(nan_obs)), ", "), behaviors{i_bhv})            
        end

        % After potential NaN removal, add single session values.
        pca_array = [pca_array; curr_bhv];        
    else
        error("Invalid value for 'observations'. Must be: 'session_mean' " + ...
              "| 'single_session'.")
    end
end

% Perform PCA.
if nv.normalize, pca_array = normalize(pca_array); end
pca_out = pca_svd(pca_array, 'corr_mat', false, 'return_all', true);

end

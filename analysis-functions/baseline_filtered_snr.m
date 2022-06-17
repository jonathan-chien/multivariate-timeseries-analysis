function [snr, mean_snr, std_snr] = baseline_filtered_snr(N, nv)
% Takes the entire baseline period and divide it into two halves; one half
% is treated as signal and the other as noise, and the SNR is optimized via
% the generalized eigenvalue problem. The baseline timeseries can be split
% in different ways repeatedly to generate a distribution of SNR.
%
% PARAMETERS
% ----------
% N -- m x n array of neural data, where m equals the number of brain
%      regions/variables and n the number of timepoints in the timeseries.
% Name-Value Pairs (nv)
%   'n_reps'        : (scalar, default = 1), specify the number of
%                     times that the baseline will be split at a random
%                     point (this is done by applying the cirshift function
%                     with a random shift and then splitting in the
%                     middle).
%   'n_dims'        : (Scalar, default = 1), dimensionality of SNR upper
%                     bound. Must not exceed number of regions/variables.
%                     See documentation in filtered_snr.m for more info.
%   'z-score'       : (1 (default) | 0), specify whether or not to zscore
%                     data before calculation of SNR upper bound.
%   'normalization' : ('bessel_corrected' | 'n_obs' | 'trace' | 'none'),
%                     specify normalization for the covariance matrix
%                     during SNR calculation.
%
% RETURNS
% -------
% snr      : p x r array, where p is the dimensionality of the calculated
%            SNR (see 'n_dims' above) and r is the number of repetitions
%            (see 'n_reps' above).
% mean_snr : p x 1 array whose i_th element is the mean over r repetitions
%            of the i_th dimension of the SNR.
% std_snr  : p x 1 array whose i_th element is the standard deviation over
%            r repetitions of the i_th dimension of the SNR.
%
% Author: Jonathan Chien 6/10/22.


% Parse inputs.
p = inputParser;
addRequired(p, 'N');
addParameter(p, 'n_reps', 1, @(x) isscalar(x) && x > 0 && floor(x) == x);
addParameter(p, 'n_dims', 1, @(x) isscalar(x) && floor(x) == x);
addParameter(p, 'z_score', true);
addParameter(p, 'normalization', 'trace');

% Preallocate.
snr = NaN(nv.n_dims, nv.n_reps);

for i_rep = 1:nv.n_reps
    % If only one repetition requested, do not perform shift.
    if nv.n_reps == 1
        [half_1, half_2] = divide_timeseries(N);
    else
        % Draw integer from {1:n_timepoints} to determine how much to shift
        % timeseries before divison for current repetition.
        shift = randi(n_timepoints, 1, 1);

        % Shift timeseries circularly. This results in a different split of
        % the entire timeseries for the SNR calculation.
        [half_1, half_2] = divide_timeseries(circshift(N), shift, 2);
    end

    % Calculate filtered SNR for half of baseline vs other half.
    snr(:,i_rep) = filtered_snr(half_1, half_2, ...
                                nv.n_dims, nv.z_score, nv.normalization);
end

% Calculate mean and std across all repetitions.
mean_snr = mean(snr, 2);
std_snr = std(snr, 0, 2);

end


% --------------------------------------------------
function [half_1, half_2] = divide_timeseries(N)

% Divide timeseries at or as near to midpoint as possible.
n_timepoints = size(N, 2);
if mod(n_timepoints, 2) == 0
    midpoint = n_timepoints/2;
else
    midpoint = ceil(n_timepoints/2);
end

half_1 = N(:,1:midpoint);
half_2 = N(:,midpoint+1:end);

end

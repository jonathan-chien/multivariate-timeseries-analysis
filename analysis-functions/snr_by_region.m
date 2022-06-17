function snr = snr_by_region(S, N, varargin)
% For an m x n signal array and m x n noise array N, where m = n_regions,
% and n = n_timepoints, calculate the signal to noise ratio separately for
% each region.
% 
% PARAMETERS
% ----------
% S : n_regions x n_timepoints signal array.
% N : n_regions x n_timepoints noise array.
% w : Optional argument; integer scalar specifying weighting scheme for
%     varaince calculation.
%
% RETURNS
% -------
% snr : n_regions x 1 vector of SNRs, where the i_th element is the SNR for
%       the i_th region.
%
% Author: Jonathan Chien 6/8/22.

% Default weighting scheme = 0, if not specified by user.
if isempty(varargin), w = 0; else, w = varargin{1}; end

snr = var(S, w, 2) ./ var(N, w, 2);

end
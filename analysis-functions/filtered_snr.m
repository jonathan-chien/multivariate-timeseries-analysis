function snr = filtered_snr(S, N, n_dims, zscore_flag, normalization, db_flag)
% Takes as input a signal array S and a signal array N and optimizes a
% filter to calculate the signal to noise ratio (SNR) upper bound. This
% filtered SNR may be multidimensional.
%
% PARAMETERS
% ----------
% S             : m x n_1 signal array, where m equals the number of
%                 variables (regions) and n_1 the number of timepoints.
% N             : m x n_2 noise array, where m equals the number of
%                 variables (regions) and n_2 the number of timepoints.
% n_dims        : Integer p subject to 0 < p < m, where m is the first dim
%                 size of S and N. p is the dimensionality of the
%                 calculated signal to noise ratio.
% zscore_flag   : (1 | 0), specify whether or not to z-score each variable
%                 (row of input arrays).
% normalization : ('bessel_corrected' | 'n_obs' | 'trace' | 'none'),
%                 specify normalization for the covariance matrix.
% db_flag       : (1 | 0), specify whether or not to convert the calculated
%                 SNR to decibels (using 10*log10).
%
% RETURNS
% -------
% snr : n_dims x 1 array whose i_th element is the scalar SNR computed
%       using the i_th eigenvector (ordered by descending eigenvalue size).
%
% Author: Jonathan Chien. 5/27/22. 


% Check inputs.
assert(size(S, 1) == size(N, 1));
assert(isscalar(n_dims));
assert(floor(n_dims) == n_dims);
assert(n_dims > 0 && n_dims <= size(S, 1));

% Option to z-score. If false, at least ensure data are mean-centered
% (row-wise) before any futher operations.
if zscore_flag
    S = zscore(S, 0, 2);
    N = zscore(N, 0, 2);
else
    S = S - mean(S, 2);
    N = N - mean(N, 2);
end

% Calculate covariance matrix. Note that if z-scoring was applied, these
% will be correlation matrices if the normalization when calculating the
% covariance matrix is by the number of degrees of freedom with the same
% correction used when z-scoring (using n - 1 for sample std here).
C_s = cov_mat(S, normalization);
C_n = cov_mat(N, normalization);

% Solve generalized eigenvalue problem. NB: eig.m does not by default
% always return eigenvalues/eigenvectors sorted by eigenvalue size. Thus
% sort eigenvalues/vectors and take top 2 dominant eigenvectors.
[V, L] = eig(C_s, C_n);
[~, ind] = sort(diag(L), 'descend');
ind = ind(1:n_dims);

% Calculate signal to noise ratio (i.e., the ratio of variances of
% optimized linear combination of signal and noise timeseries,
% respectively).
snr = diag((V(:,ind)' * C_s * V(:,ind)) ./ (V(:,ind)' * C_n * V(:,ind)));

% Optionally convert SNR to dB units.
if db_flag, snr = 10 * log10(snr); end

end


% --------------------------------------------------
function C = cov_mat(X, normalization)
% Calculate n x n covaraince matrix for n vectors of length t passed in as
% an n x t array. Arg normalization specifies how to normalize the
% covariance matrix.

% Redundant here but better to include.
X_m = X - mean(X, 2);

% Calculate covariance matrix.
C = X_m * X_m';

% Apply specified normalization.
switch normalization
    case 'bessel_corrected'
        norm = size(X, 2) - 1;
    case 'n_obs'
        norm = size(X, 2);
    case 'trace'
        norm = trace(C);
    case 'none'
        norm = 1;
    otherwise
        error('Please provide a valid normalization option.')
end
 
C = C / norm;

end

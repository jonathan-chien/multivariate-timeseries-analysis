function null = get_null_timeseries(timeseries, method, varargin)
% Accepts as input an n_vars x n_timepoints array and generates a null
% distribution for each variable/the set of variables.
%
% PARAMETERS
% ----------
% timeseries    : n_vars x n_timepoints array.
% method        : ('permutation' (default) | 'jitter'), specify which
%                 method used to generate null draw. The default
%                 option, 'permutation', permutes the timepoints for
%                 each variable independently, suitable for general
%                 analyses. The 'jitter' method leaves each variable's
%                 timeseries largely intact but applies a circular
%                 shift, suitable for analyses such as
%                 cross-correlation that take into account the lag
%                 between timeseries, or for covariance/correlation
%                 based analyses in general.
% jitter_offset : (Optional), If 'method' = 'jitter', pass a scalar which
%                 controls the strength of the jittering. For
%                 'jitter_offset' = 20 (default), the jitter for each
%                 variable will be drawn uniformly from the set of integers
%                 falling on the inclusive interval [-20 20]. This argument
%                 is ignored if method = 'permutation'.
%
% RETURNS
% -------
% null : n_vars x n_timepoints null array.
%
% Author: Jonathan Chien. Refactored from body of xcorr_by_bhv on 5/31/22.


% Check if user passed optional jitter offset argument. Set to 10 if not.
jitter_offset = varargin{1};
if isempty(jitter_offset), jitter_offset = 10; end

% Array sizes/preallocation.
[n_vars, n_timepts] = size(timeseries);
null = NaN(size(timeseries));

% Generate one draw from underlying null.
switch method
    case 'permutation'
        % Permute each region's timeseries independently.    
        for i_var = 1:n_vars
            null(i_var,:) = timeseries(i_var,randperm(n_timepts));
        end
    case 'jitter'
        % Jitter each region's timeseries independently.
        for i_var = 1:n_vars
            % Flip current variable's jitter to negative displacement with
            % prob 0.5.
            draw = randi(2, 1, 1); 
            if mod(draw, 2) == 0, flip = -1; else, flip = 1; end

            % Apply jitter.
            null(i_var,:) = circshift(timeseries(i_var,:), ...
                                      randi(jitter_offset, 1, 1) * flip, ...
                                      2);
        end
    otherwise
        error("Specify method as 'permutation' or 'jitter'.")
end

end

function out = xcorr_by_bhv(timeseries, session_raw, behavior, varargin)
% Accepts a multivariate timeseries and extracts portions
% corrresponding to behavior of interest and calculates cross-correaltions
% at different lags, with significance measures available via permutation
% or jittering. Function also provides options for visualization.
%
% PARAMETERS
% ----------
% timeseries  : n_regions x n_timepoints array of neural activity. This is
%               passed separately (rather than using, e.g., 
%               session_raw.Lfold) in case we want to analyze deconvolved
%               data instead of raw calcium/photometry traces.
% session_raw : Scalar struct in one cell out of the cell array output of
%               extract_single_subj). This struct contains a list of all
%               behaviors form the session (needed to check that the
%               specified behavior occured during the passed in session)
%               and the start and stop times of all epochs.
% behavior    : String value specifiying the behavior of interest. Note
%               capitalization.
% Name-Value Pairs (nv) 
%   'bin_width'      : (Scalar) number of timepoints in a bin.
%   'step_size'      : (Scalar) step size by which each bin will be
%                      advanced. Units are in timepoints.
%   'zscore'         : (1 (default) | 0), specify whether or not to zscore
%                      data after extracting the behavior of interest and
%                      baseline.
%   'scaleopt'       : (String: 'none' (default) | 'biased' | 'unbiased' |
%                      'normalized' | 'coeff'), specify normalization
%                      option for xcorr.m.
%   'maxlag'         : (Scalar), specify the maximum number of lags to be
%                      computed using xcorr.m
%   'svn_opts'       : 2-vector whose first and second elements correspond
%                      respectively to the 'normalize' and 'scale' options
%                      for the local function calc_svn. If 'normalize' is
%                      true, the difference between signal and noise will
%                      be normalized by the noise. If 'scale' is true, the
%                      output of xcorr.m (across correlation pairs and
%                      lags) will be scaled to have standard deviation 1.
%   'pval'           : ('right-tailed' (default) | 'left-tailed' |
%                      'two-tailed' | false). If false, calculation of null
%                      distribution will be suppressed. If not false,
%                      specify the type of test desired.
%   'null_method'    : ('permutation' (default) | 'jitter'), specify
%                      the method used to generate each draw from the null
%                      distribution. If 'permutation' the timeseries for
%                      each brain region (each row of input argument
%                      "timeseries") will be permuted indepednently. If
%                      'jitter', each brain region's timeseries will be
%                      indpendently shifted using a random shift passed to
%                      the circshift function.
%   'jitter_offset'  : (Scalar, default = 10), if 'null_method' =
%                      'jitter', this argument controls the range of values
%                      from which we randomly choose when determining how
%                      many much (i.e., by how many timepoints) to jitter a
%                      region's timeseries. If x is the value passed in,
%                      values are chosen randomly from the set of integers
%                      {1, 2, ... , x}. A different random draw is made for
%                      each region, for each draw from the null. If
%                      'null_method' does not equal 'jitter', this argument
%                      is ignored.
%   'n_null'         : (Scalar, default = 1000), specify the number
%                      of null datasets. Each dataset is generated by one
%                      application of the permutation or jitter methods
%                      described above; all calculations are performed
%                      on each null datset, resulting in a null
%                      distribution for each metric, with the number
%                      of draws equaling 'n_null'.
%   'return_null'    : (1 | 0 (default), specify whether or not to return
%                      the null distribution for each metric.
%   'extract_opt_by' : ('smallest_p_val' | 'max_xcorr' (default)), specify
%                      how to extract optimal lags. If 'max_xcorr', the lag
%                      resulting in the highest cross-correlation will be
%                      chosen for each cross-correlation pair is chosen. If
%                      'smallest_p_val', the lag resulting in the smallest
%                      p value will be chosen for each cross-correlation
%                      pair.
%   'plot_all_lags'  : (1 (default) | 0), for each lag, extract
%                      non-redundant elements of each correlation matrix
%                      (lower triangular including autocorrelations on main
%                      diagonal) and flatten; then the plot matrix where
%                      the i_th row has the flattened extracted elements
%                      for the i_th lag. R (behavior), R_0 (baseline), SVN
%                      (signal vs noise), as well as p values for these
%                      three are all plotted.
%   'plot_opt_lag'   : ('flattened' (default) | 'corr_mat' | 'adjacency'),
%                      specify the format to plot the cross-correlations at
%                      the optimal chosen lag (see 'extract_opt_by'). If
%                      'flattened', the flattened correlation matrix (with
%                      redundant correlation pairs removed) will be
%                      plotted. If 'corr_mat', a cross-correlation matrix
%                      will be plotted. If 'adjacency', an adjacency matrix
%                      will be plotted with a 1 where the cross-correlation
%                      between regions is significant and a 0 otherwise.
%   'threshold'      : (false (default) | scalar), specify p value
%                      threshold to be used to identify significant
%                      regional interactions and calculate adjacencey
%                      matrix.
%   'cmap'           : (colormap name, default = parula)
%
% RETURNS
% -------
% out : Scalar struct with the following fields (where m = 2 * max_lag + 1
%       and n = square of the number of regions):
%   .R         : m x n array whose i_th, j_th element is the correlation of
%                the j_th pair (equivalent to linear index of corresponding
%                element of a correlation matrix) at the i_th lag during
%                the behavior of interest. If a p value threshold was set
%                (see 'threshold' under PARAMETERS), all elements whose p
%                value was above the threshold will be set to NaN.
%   .R_0       : m x n array whose i_th, j_th element is the correlation of
%                the j_th pair (equivalent to linear index of corresponding
%                element of a correlation matrix) at the i_th lag during
%                the baseline period. If a p value threshold was set (see
%                'threshold' under PARAMETERS), all elements whose p value
%                was above the threshold will be set to NaN.
%   .svn       : m x n array whose i_th, j_th element is the comparison in
%                correlation of the j_th pair (equivalent to linear index
%                of corresponding element of a correlation matrix) at the
%                i_th lag between the behavior and baseline period (see
%                'svn_opts' under PARAMETERS). If a p value threshold was
%                set (see 'threshold' under PARAMETERS), all elements whose
%                p value was above the threshold will be set to NaN.
%   .p         : Scalar struct with fields (not returned if 'pval' =
%                false):
%       .R   : m x n array whose i_th, j_th element is the p value
%              attaching to the i_th, j_th element of out.R.
%       .R_0 : m x n array whose i_th, j_th element is the p value
%              attaching to the i_th, j_th element of out.R_0.
%       .svn : m x n array whose i_th, j_th element is the p value
%              attaching to the i_th, j_th element of out.svn.
%   .null      : Scalar struct with the following fields:
%       .R   : m x n x p array whose i_th, j_th, k_th element is the k_th
%              draw from the null distribution for the i_th, j_th element
%              of out.R.
%       .R_0 : m x n x p array whose i_th, j_th, k_th element is the k_th
%              draw from the null distribution for the i_th, j_th element
%              of out.R_0.
%       .svn : m x n x p array whose i_th, j_th, k_th element is the k_th
%              draw from the null distribution for the i_th, j_th element
%              of out.svn.
%   .adjacency : Scalar struct with the following fields (not returned if
%                'pval' = false):
%       .R   : m x n array whose i_th, j_th element has a value of 1 if the
%              i_th, j_th element of out.p.R was below the value set by
%              'threshold' (see PARAMETERS) and a 0 otherwise.
%       .R_0 : m x n array whose i_th, j_th element has a value of 1 if the
%              i_th, j_th element of out.p.R_0 was below the value set by
%              'threshold' (see PARAMETERS) and a 0 otherwise.
%       .svn : m x n array whose i_th, j_th element has a value of 1 if the
%              i_th, j_th element of out.p.svn was below the value set by
%              'threshold' (see PARAMETERS) and a 0 otherwise.
%   .opt       : Scalar struct with the following fields:
%       .R         : 1 x n array whose j_th element is the behavioral
%                    cross-correlation at the optimal lag, for the j_th
%                    region pair (optimal value of the j_th column of
%                    out.R, as defined by the 'extract_opt_by' name-value
%                    pair (see PARAMETERS)).
%       .R_0 :     : 1 x n array whose j_th element is the baseline
%                    cross-correlation at the optimal lag, for the j_th
%                    region pair (optimal value of the j_th column of
%                    out.R_0, as defined by the 'extract_opt_by' name-value
%                    pair (see PARAMETERS)).
%       .svn       : 1 x n array whose j_th element is the comparison of
%                    behavioral vs baseline cross-correlation at the
%                    optimal lag, for the j_th region pair (optimal value
%                    of the j_th column of out.svn, as defined by the
%                    'extract_opt_by' name-value pair (see PARAMETERS)).
%       .p         : Scalar struct with the following fields (not returned
%                    if 'pval' = false):
%           .R   : 1 x n array whose j_th element is the p value for the
%                  behavioral cross-correlation at the optimal (as defined
%                  by the 'extract_opt_by' name-value pair (see
%                  PARAMETERS)) lag, for the j_th region pair.
%           .R_0 : 1 x n array whose j_th element is the p value for the
%                  baseline cross-correlation at the optimal (as defined
%                  by the 'extract_opt_by' name-value pair (see
%                  PARAMETERS)) lag, for the j_th region pair.
%           .svn : 1 x n array whose j_th element is the p value for the
%                  comparison of signal vs noise cross-correlation at the
%                  optimal (as defined by the 'extract_opt_by' name-value
%                  pair (see PARAMETERS)) lag, for the j_th region pair.
%                  (see as well, 'svn_opts' under PARAMETERS).
%       .adjacency : Scalar struct with the following fields (not returned
%                    if 'pval' = false):
%           .R   : 1 x n array whose j_th element has a value of 1 if the p
%                  value for the behavioral cross-correlation at the
%                  optimal (as defined by the 'extract_opt_by' name-value
%                  pair (see PARAMETERS)) lag, for the j_th region pair,
%                  was below the p value threshold set via the 'threshold'
%                  name-value pair, and a 0 otherwise.
%           .R_0 : 1 x n array whose j_th element has a value of 1 if the p
%                  value for the baseline cross-correlation at the
%                  optimal (as defined by the 'extract_opt_by' name-value
%                  pair (see PARAMETERS)) lag, for the j_th region pair,
%                  was below the p value threshold set via the 'threshold'
%                  name-value pair, and a 0 otherwise.
%           .svn : 1 x n array whose j_th element has a value of 1 if the p
%                  value for the comparison of behavioral vs baseline
%                  cross-correlation at the optimal (as defined by the
%                  'extract_opt_by' name-value pair (see PARAMETERS)) lag,
%                  for the j_th region pair, was below the p value
%                  threshold set via the 'threshold' name-value pair, and a
%                  0 otherwise.
%
% Author: Jonathan Chien. 5/20/22. Version 0.2.


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addRequired(p, 'behavior', @ischar);
addParameter(p, 'bin_width', 5, @(x) isscalar(x) || floor(x) == x); % This is in units of timepoints
addParameter(p, 'step_size', 2, @(x) isscalar(x) || floor(x) == x)
addParameter(p, 'zscore', true, @islogical); % Z-score data before calling xcorr
addParameter(p, 'scaleopt', ...
             @(x) ismember(x, {'none', 'biased', 'unbiased', 'normalized', 'coeff'}))
addParameter(p, 'maxlag', 10, @(x) isscalar(x) || floor(x) == x);
addParameter(p, 'svn_opts', [true true]); % Vector whose elements are values for the nvp of calc_svn
addParameter(p, 'pval', 'right-tailed');
addParameter(p, 'null_method', 'permutation', @(x) ismember(x, {'permutation', 'jitter'}))
addParameter(p, 'jitter_offset', 10, @isscalar);
addParameter(p, 'n_null', 1000, @isscalar);
addParameter(p, 'return_null', false, @islogical);
addParameter(p, 'extract_opt_by', 'max_xcorr');
addParameter(p, 'plot_all_lags', true, @islogical);
addParameter(p, 'plot_opt_lag', 'flattened');
addParameter(p, 'threshold', false); % Either false, or a p value threshold
addParameter(p, 'cmap', parula);
parse(p, timeseries, session_raw, behavior, varargin{:});
nv = p.Results; 


%% Preprocessing

% Each region consists of a scalar timeseries.
n_regions = size(timeseries, 1);

% Pull out parts of timeseries labeled as belonging to behavior of
% interest.
timeseries_bhv = isolate_behavior(timeseries, session_raw, behavior);

% Pull out parts of timeseries labeled as baseline.
timeseries_base = isolate_baseline(timeseries, session_raw);

% Bin data.
bin_timeseries_bhv = bin_data(timeseries_bhv, nv.bin_width, nv.step_size);
bin_timeseries_base = bin_data(timeseries_base, nv.bin_width, nv.step_size);

% Option to z-score data first.
if nv.zscore
    bin_timeseries_bhv = normalize(bin_timeseries_bhv, 2);
    bin_timeseries_base = normalize(bin_timeseries_base, 2);
end


%% Cross-correlation 

% Calculate cross-correlation matrix for targeted and baseline behavior.
[out.R, lags] = xcorr(transpose(bin_timeseries_bhv), nv.maxlag, nv.scaleopt);
out.R_0 = xcorr(transpose(bin_timeseries_base), nv.maxlag, nv.scaleopt);

% Calculate signal to noise metrics on actual data.
out.svn = calc_svn(out.R, out.R_0, ...
                   'normalize', nv.svn_opts(1), 'scale', nv.svn_opts(2));

% Get names of fields of out before further additions.
mat_names = fieldnames(out);

% Number of lags used.
n_lag = length(lags);


%% Generate null distribution 

if nv.pval

% Preallocate.
n_corr = n_regions^2;
for i_field = 1:length(mat_names)
    null.(mat_names{i_field}) = NaN(n_lag, n_corr, nv.n_null);
end

% Generate null distribution using permuted/jittered data.
for i_null = 1:nv.n_null
    % Generate one draw from underlying null.
    null_bin_timeseries_bhv ...
        = get_null_timeseries(bin_timeseries_bhv, nv.null_method, nv.jitter_offset);
    null_bin_timeseries_base ...
        = get_null_timeseries(bin_timeseries_base, nv.null_method, nv.jitter_offset);
    
    % Calculate cross-correlation on null.
    null.R(:,:,i_null) ...
        = xcorr(transpose(null_bin_timeseries_bhv), nv.maxlag, nv.scaleopt);
    null.R_0(:,:,i_null) ...
        = xcorr(transpose(null_bin_timeseries_base), nv.maxlag, nv.scaleopt);

    % Calculate SVN metric.
    null.svn(:,:,i_null) = calc_svn(null.R(:,:,i_null), null.R_0(:,:,i_null), ...
                                    'normalize', nv.svn_opts(1), ...
                                    'scale', nv.svn_opts(2));
end

% Preallocate.
for i_field = 1:length(mat_names)
    out.p.(mat_names{i_field}) = NaN(n_lag, n_regions^2);
end

% Attach p values. 
for i_lag = 1:2*nv.maxlag + 1
for i_field = 1:length(mat_names)
    out.p.(mat_names{i_field})(i_lag,:) ...
        = tail_prob(out.(mat_names{i_field})(i_lag,:)', ...
                    squeeze(null.(mat_names{i_field})(i_lag,:,:)), ...
                    2, ...
                    'type', nv.pval, 'exact', false, 'correction', true);
end
end

% Option to return null distributions.
if nv.return_null, out.null = null; end

end


%% Optionally apply a threshold and get adjacency matrix

% Everything below threshold is set to 0. Adjacency matrix is created by
% assigning 1 if threshold is passed and 0 otherwise.
if nv.threshold & nv.pval
    for i_field = 1:length(mat_names)
        out.(mat_names{i_field})(out.p.(mat_names{i_field}) >= nv.threshold) = NaN;
        out.adjacency.(mat_names{i_field}) ...
            = out.p.(mat_names{i_field}) < nv.threshold;
    end
end


%% Extract max cross-correlation (i.e., find optimal lag)

% For each interaction, take the optimal lag, either by the largest
% cross-correlation across lags, or the smallest p value across lags (if p
% values were computed).
for i_field = 1:length(mat_names)
    switch nv.extract_opt_by
        case 'max_xcorr'
            % Extract maximum values (R, R_0, SVN). 
            [out.opt.(mat_names{i_field}), opt_ind] ...
                = max(out.(mat_names{i_field}), [], 1, 'ComparisonMethod','abs');
        case 'smallest_p_val'
            % Throw error if p values were not computed.
            assert(any(nv.pval), ...
                   "'extract_opt_by' may not be 'smallest_p_val' if " + ...
                   "p values were not computed ('pval' = false).")
            [out.opt.(mat_names{i_field}), opt_ind] ...
                = min(out.p.(mat_names{i_field}), [], 1);
    end
    
    % Extact corresponding p values (if calculated) and adjacency values as
    % well.
    if nv.pval
        out.opt.p.(mat_names{i_field}) = NaN(1, n_regions^2);
        out.opt.adjacency.(mat_names{i_field}) = NaN(1, n_regions^2);
        for i_col = 1:n_regions^2
            out.opt.p.(mat_names{i_field})(i_col) ...
                = out.p.(mat_names{i_field})(opt_ind(i_col),i_col);
            out.opt.adjacency.(mat_names{i_field})(i_col) ...
                = out.adjacency.(mat_names{i_field})(opt_ind(i_col),i_col);
        end
    end
end


%% Optional plotting

% Plot flattened correlation matrices, but only the lower triangular
% portion of the matrix (main diagonal optional). See documentation for
% more information.
if nv.plot_all_lags

    % Set whether or not to include autocorrelations (on main diagonal of
    % correlation matrix).
    MAIN_DIAG = false;
    
    % Get subscripts of correlation matrices.
    subscripts = get_subscripts(n_regions, MAIN_DIAG);
    
    % Plot output of xcorr with redundancies removed.
    figure   
    h = cell(length(mat_names), 1);
    for i_field = 1:length(mat_names)
        % Create heatmap for current matrix and set properties.
        subplot(length(mat_names), 1, i_field)
        h{i_field} = heatmap( ...
            extract_lower_flattened(out.(mat_names{i_field}), MAIN_DIAG) ...
                             );       
        config_all_lags_flattened_heatmap(h{i_field}, nv.cmap, subscripts, lags);

        title(set_title(mat_names{i_field}, behavior))
    end

    sgtitle('Cross-correlations')


    % If p values were computed, plot them. 
    if nv.pval
        figure
        h = cell(length(mat_names), 1);

        for i_field = 1:length(mat_names)
            % Create heatmap for current matrix and set properties .
            subplot(3, 1, i_field)
            h{i_field} = heatmap( ...
                extract_lower_flattened(out.p.(mat_names{i_field}), MAIN_DIAG) ...
                                 );
            config_all_lags_flattened_heatmap(h{i_field}, nv.cmap, subscripts, lags);
    
            % Set subplot title.
            title(set_title(mat_names{i_field}, behavior))
        end

        sgtitle('P values')
    end

end
    

% For each interaction between areas (no redundancy), take max value across
% lags. Plot each flattened matrix as a row.
if strcmp(nv.plot_opt_lag, 'flattened')

    % Set whether or not to include autocorrelations (on main diagonal of
    % correlation matrix).
    MAIN_DIAG = false;
    
    % Get subscripts of correlation matrices.
    subscripts = get_subscripts(n_regions, MAIN_DIAG);

    % Get max cross-correlation vals and their corresponding p values.
    % Optionally sort interactions by the largest difference between
    % baseline and behavior of interest (regardless of whether
    % normalization by baseline was applied).
    opt_xcorr_vals = [extract_lower_flattened(out.opt.R, MAIN_DIAG); ...
                      extract_lower_flattened(out.opt.R_0, MAIN_DIAG); ...
                      extract_lower_flattened(out.opt.svn, MAIN_DIAG)];
    opt_p_vals = [extract_lower_flattened(out.opt.p.R, MAIN_DIAG); ...
                  extract_lower_flattened(out.opt.p.R_0, MAIN_DIAG); ...
                  extract_lower_flattened(out.opt.p.svn, MAIN_DIAG)];
    if true
        [~,diff_ind] = sort(abs(diff(opt_xcorr_vals(1:2,:))), 'descend');
        opt_xcorr_vals = opt_xcorr_vals(:,diff_ind);
        opt_p_vals = opt_p_vals(:,diff_ind);
        subscripts = subscripts(diff_ind);
    end

    % Plot cross-correlations.
    figure
    h = heatmap(opt_xcorr_vals);
    config_opt_lag_flattened_heatmap(h, nv.cmap, subscripts, behavior);
    title('Cross-correlations')
  
    % If p values were computed, plot them.
    if nv.pval
        figure
        h = heatmap(opt_p_vals);
        config_opt_lag_flattened_heatmap(h, nv.cmap, subscripts, behavior);     
        title('P values')
    end

% For each interaction between areas (no redundancy), take max value across
% lags, then reshape into matrix format.
elseif any(strcmp(nv.plot_opt_lag, {'corr_mat', 'adjacency'}))

    % Set indices to be excluded using FULL_MAT. Options include either
    % lower triangular portion only (if FULL_MAT = false) or full matrix
    % (if FULL_MAT = true). In either case, use MAIN_DIAG to specify
    % whether or not to include autocorrelations (on main diagonal of
    % matrix).
    MAIN_DIAG = false;
    FULL_MAT = false;
    if ~FULL_MAT % Show only lower triangular portion
        exclusion_ind ...
            = setdiff(1:n_regions^2, get_lower_tri_linear_ind(n_regions, MAIN_DIAG));
    else % Show full matrix.
        if MAIN_DIAG % Include main diagonal.
            exclusion_ind = []; 
        else % Exclude main diagonal.
            exclusion_ind = 1:n_regions+1:n_regions^2;
        end
    end

    % Reshape from flattened into matrix form. Set non lower triangular
    % values to NaN.
    figure
    h = cell(length(mat_names), 1);
    for i_field = 1:length(mat_names)
        subplot(length(mat_names), 1, i_field)

        % Reshape from flattened.
        switch nv.plot_opt_lag
            case 'corr_mat'
                mat = reshape(out.opt.(mat_names{i_field}), [n_regions n_regions]);
            case 'adjacency'
                mat = reshape(out.opt.adjacency.(mat_names{i_field}), [n_regions n_regions]);
        end
        
        % Disqualify elements as described above and create heatmap.
        mat(exclusion_ind) = NaN; 
        
        % Create heatmap and set heatmap properties.
        h{i_field} = heatmap(mat);
        config_opt_lag_matrix_heatmap(h{i_field}, nv.cmap, 'none', nv.plot_opt_lag);
        title(set_title(mat_names{i_field}, behavior))
    end

    if strcmp(nv.plot_opt_lag, 'adjacency')
        sgtitle('Cross-correlations adjacency matrices')
    else
        sgtitle('Cross-correlations')
    end
    
    % Plot p values if they were computed. Note that there are no p values
    % associated with the adjacency matrices, so this option is valid only
    % for the 'plot_opt_lag', 'corr_mat' name-value pair.
    if nv.pval & strcmp(nv.plot_opt_lag, 'corr_mat')

        figure
        h = cell(length(mat_names), 1);
        
        for i_field = 1:length(mat_names)
            subplot(length(mat_names), 1, i_field)

            % Reshape from flattened, disqualify elements as above, and
            % create heatmap.
            mat = reshape(out.opt.p.(mat_names{i_field}), [n_regions n_regions]);
            mat(exclusion_ind) = NaN;
            h{i_field} = heatmap(mat);
            
            % Set heatmap properties.
            config_opt_lag_matrix_heatmap(h{i_field}, nv.cmap, 'none', nv.plot_opt_lag);
            
            title(set_title(mat_names{i_field}, behavior))
        end
        
        sgtitle('P values')
    end

end

end


%% Local functions

% --------------------------------------------------
function svn = calc_svn(R, R_0, varargin)
% Calculate signal to noise measure comparing cross-correlation of regions
% during periods of interest compared with baseline cross-correlation.

p = inputParser;
addRequired(p, 'R');
addRequired(p, 'R_0');
addParameter(p, 'normalize', true, @islogical);
addParameter(p, 'scale', true, @islogical); % Divide svn by total std across all elements; this is done as the final step
parse(p, R, R_0, varargin{:})

% For each interaction at each lag, calculate difference.
svn = R - R_0;

% Optionally normalize (for each interaction at each lag).
if p.Results.normalize, svn = svn ./ R_0; end

% Optionally scale final matrix by standard deviation.
if p.Results.scale
    stdev = std(svn, 0, 'all');
    svn = svn / stdev;
end

end

% --------------------------------------------------
function tril_extracted = extract_lower_flattened(R, main_diag_flag)
% Input R is an m x n^2 array where the i_th row contains the flattened n x
% n correlation matrix with the i_th lag. main_diag_flag is passed to local
% function get_lower_tri_linear_ind (see it for more info).

% Get array sizes and linear indices of elements on and below main diagonal
% of correlation matrix.
n_lag = size(R, 1);
n_vars = sqrt(size(R, 2));
if main_diag_flag
    n_extracted = n_vars * (n_vars + 1) / 2;
elseif ~main_diag_flag
    n_extracted = n_vars * (n_vars - 1) / 2;
end
lower_tri_linear_ind = get_lower_tri_linear_ind(n_vars, main_diag_flag);

% Extract elements.
tril_extracted = NaN(n_lag, n_extracted);
for i_lag = 1:n_lag
    tril_extracted(i_lag,:) = R(i_lag,lower_tri_linear_ind);
end

end


% --------------------------------------------------
function linear_ind = get_lower_tri_linear_ind(n, main_diag_flag)
% Get linear indices of all elements on and below main diagonal of an n x n
% matrix. If main_diag_flag == 1, the indices of the main digonal elements
% will be included. If main_diag_flag = 0, indices of the main digonals are
% excluded.

i_col = 1;
linear_ind = [];
while i_col <= n
    to_add = (1:n) + ((i_col - 1) * n);
    if main_diag_flag % Include main diagonal
        to_add = to_add(i_col:end);
    elseif ~main_diag_flag % Exclude main diagonal
        to_add = to_add(i_col+1:end);
    end
    linear_ind = [linear_ind, to_add];
    i_col = i_col + 1;
end

end


% --------------------------------------------------
function subscripts = get_subscripts(n, main_diag_flag)
% Get lower triangular subscripts of n x n matrix (with main
% diagonal if main_diag_flag == 1 and without if main_diag_flag == 0).

[row, col] = ind2sub([n n], get_lower_tri_linear_ind(n, main_diag_flag));

subscripts = cell(length(row), 1);
for i_idx = 1:length(row)
    subscripts{i_idx} = strjoin(string([row(i_idx) col(i_idx)]), ',');
end

end


% --------------------------------------------------
function config_all_lags_flattened_heatmap(h, cmap, subscripts, lags)
% Configure heatmap for flattened correlation matrices indexed by lag.

% If true, set the angle of the X tick labels manually. Otherwise they are
% set dynamically based on figure size etc.
SET_ANGLE = false;

h.GridVisible = false;
h.Colormap = cmap;
h.ColorbarVisible = true;
h.MissingDataLabel = 'Below threshold';

% Set X labels.
h.XDisplayLabels = subscripts;

% Prepare and set Y labels.
yt = num2cell(lags);
yt(setdiff(1:length(lags), floor(linspace(1, length(lags), 5)))) = {''};
h.YDisplayLabels = yt;

% Option to manually set X tick angles.
if SET_ANGLE
    s = struct(h);
    s.XAxis.TickLabelRotation = 45;
end

end


% --------------------------------------------------
function config_opt_lag_flattened_heatmap(h, cmap, subscripts, behavior)
% Configure heatmap for flattened correlation matrices with elements
% consisting of maximum values across all lags (opt is for optimal).

% If true, set the angle of the X tick labels manually. Otherwise they are
% set dynamically based on figure size etc.
SET_ANGLE = false;

h.GridVisible = false;
h.Colormap = cmap;
h.ColorbarVisible = true;
h.MissingDataLabel = 'Below threshold';

% Set X labels.
h.XDisplayLabels = subscripts;

% Prepare and set Y labels.
h.YDisplayLabels = {behavior, 'baseline', 'Signal vs noise'};

% Option to manually set X tick angles.
if SET_ANGLE
    s = struct(h);
    s.XAxis.TickLabelRotation = 45;
end

end


% --------------------------------------------------
function config_opt_lag_matrix_heatmap(h, cmap, text_flag, mat_type)
% Configure heatmap for non-flattened (square) correlation matrices at
% maximum values across all lags (that is at the optimal lag; opt =
% optimal). The mat_type argument enacts different settings for the
% correlation matrix ('corr_mat') vs the adjacency matrix 'adjacency').

h.GridVisible = false;
h.CellLabelColor = text_flag; % Controls whether or not data are printed on each cell of heatmap
h.MissingDataLabel = '';

switch mat_type
    case 'corr_mat'
        h.Colormap = cmap;
        h.ColorbarVisible = true;
        h.MissingDataColor = ones(1, 3);
    case 'adjacency'
        h.Colormap = gray;
        h.ColorbarVisible = false;
        h.MissingDataColor = ones(1, 3) * 0.5;
end

end


% --------------------------------------------------
function title = set_title(mat_name, behavior)
% Set readable titles based on matrix.

switch mat_name
    case 'svn'
        title = 'Signal vs noise';
    case 'R'
        title = sprintf('%s epochs', behavior);
    case 'R_0'
        title = 'Baseline epochs';
end

end

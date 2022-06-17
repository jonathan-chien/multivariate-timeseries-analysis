function snr = snr_by_bhv(timeseries, session_raw, behavior, varargin)
% This function accepts as input a multivariate timeseries (consisting of a
% set of scalar timeseries, one from each brain region) and a list of
% behavioral variables to test. Parts of the behavior to be tested are
% pulled out and compared against the baseline period to compute various
% kinds of signal to noise ratio. The user can bin the data by passing in
% the desired bin and step size.
%
% PARAMETERS
% ----------
% timeseries  : m x n array of neural data (preprocessing should be done
%               if desired), where m is the number of regions/variables,
%               and n the number of timepoints in the timeseries.
% session_raw : One cell of the sessions_raw variable returned by
%               extract_single_subj. Scalar struct with fields: 
%   .Lfold     : Recorded GCaMP traces
%   .t         : Time stamp of each GCaMP imaging frame
%   .behaviors : Annotated behavior
%   .Fstart    : Starting frame of behaviors
%   .FStop     : Stopping frame of behaviors
%   .regions   : Brain region
%   .FL        : Sampling rate, 25 Hz --> 40 ms
% behavior : String value that is the name of the behavior to be isolated.
% Name-Value Pairs (nv)
%   'min_timeseries_length' : (Scalar) minimum number of timepoints
%                             required for specified behavior. If the
%                             session had fewer timepoints for this
%                             behavior than the minimum number, an
%                             exception will be thrown.
%   'bin_width'             : (Scalar) number of timepoints in a bin.
%   'step_size'             : (Scalar) step size by which each bin will be
%                             advanced. Units are in timepoints.
%   'snr_dim'               : (Scalar) Number of components of the
%                             calculated SNR upper bound.
%   'zscore'                : (1 (default) | 0), specify whether or not to
%                             z-score the extracted behavior and baseline
%                             timeseries prior to SNR calculations.
%   'normalization'         : ('bessel_corrected' | 'n_obs' | 'trace' |
%                             'none'), specify normalization for the
%                             covariance matrix during SNR optimization.
%   'decibels'              : (1 (default) | 0), specify whether or not to
%                             convert raw SNR to decibels using dB =
%                             10*log10(SNR).
%   'pval'                  : ('right-tailed' (default) | 'left-tailed' |
%                             'two-tailed' | false). If false, calculation
%                             of null distribution will be suppressed. If
%                             not false, specify the type of test desired.
%   'null_method'           : ('permutation' (default) | 'jitter'), specify
%                             the method used to generate each draw from
%                             the null distribution. If 'permutation' the
%                             timeseries for each brain region (each row of
%                             input argument "timeseries") will be permuted
%                             indepednently. If 'jitter', each brain
%                             region's timeseries will be indpendently
%                             shifted using a random shift passed to the
%                             circshift function.
%   'jitter_offset'         : (Scalar, default = 10), if 'null_method' =
%                             'jitter', this argument controls the range of
%                             values from which we randomly choose when
%                             determining how many much (i.e., by how many
%                             timepoints) to jitter a region's timeseries.
%                             If x is the value passed in, values are
%                             chosen randomly from the set of integers {1,
%                             2, ... , x}. A different random draw is made
%                             for each region, for each draw from the null.
%                             If 'null_method' does not equal 'jitter',
%                             this argument is ignored.
%   'n_null'                : (Scalar, default = 1000), specify the number
%                             of null datasets. Each dataset is generated
%                             by one application of the permutation or
%                             jitter methods described above; all SNR
%                             calculations are performed on each null
%                             datset, resulting in a null distribution for
%                             each SNR measure, with the number of draws
%                             equaling 'n_null'.
%   'null_int_size'         : (Scalar, default = 95), specify the size of
%                             the interval to be computed on each null
%                             distribution.
%   'return_null'           : (1 | 0 (default), specify whether or not to
%                             return the null distribution for each SNR
%                             measure.
%
% RETURNS
% -------
% snr : Scalar struct with the following fields:
%   .by_region                 : n_regions x 1 vector whose i_th element is
%                                the SNR (var(signal) / var(noise)) for the
%                                i_th region/variable.
%   .full                      : p x 1 vector whose i_th element is the
%                                i_th component of the SNR upper bound. p
%                                is the dimensionality of the SNR upper
%                                bound.                    
%   .leave_one_out             : n_regions x p array whose i_th, j_th
%                                element is the j_th component of the SNR
%                                upper bound calculated by leaving out the
%                                i_th region.
%   .leave_one_out_diff        : n_regions x p array whose i_th, j_th
%                                element is the difference between the j_th
%                                component of the SNR upper bound
%                                calculated on all regions and the j_th
%                                component of the SNR upper bound
%                                calculated by leaving out the i_th region.
%   .leave_one_out_diff_normed : n_regions x p array whose i_th, j_th
%                                is the i_th, j_th element of
%                                .leave_one_out_diff, normed by the j_th
%                                component of the SNR upper bound using all
%                                regions (i.e., .full).
%   .sig                       : This field is present only if the 'p_val'
%                                is not false and has a valid value. If
%                                valid, sig is a scalar struct with the
%                                following fields:
%       .p         : Scalar struct of p values whose fieldnames match those
%                    of snr (listed above), except for 'sig'.
%                    Shapes and elements of arrays in each field
%                    also match and correspond one to one; e.g.,
%                    the i_th, j_th element of snr.sig.pdistributions is the p
%                    value associated with the i_th, j_th element of
%                    snr.leave_one_out.
%       .null_int  : Scalar struct of null intervals whose fieldnames match
%                    those of snr (listed above), except for 'sig'. For
%                    each field of snr, the corresponding field in
%                    sig.null_int has an array with one extra dimension (to
%                    store the lower and upper interval bounds). For
%                    example, if snr.leave_one_out has shape 13 x 2, then
%                    snr.sig.null_int_leave_one_out has shape 13 x 2 x 2,
%                    where the (i,j,:) slice contains the upper and lower
%                    interval bounds for the i_th, j_th value of
%                    snr.leave_one_out.
%       .z         : Scalar struct of z scores whose fieldnames match
%                    those of snr (listed above), except for 'sig'. Shapes
%                    and elements of arrays in each field also match and
%                    correspond one to one; e.g., the i_th, j_th element of
%                    snr.sig.z.by_region is the z-score on the null
%                    permutation distribution associated with the i_th,
%                    j_th element of snr.by_region.
%       .null_mean : Scalar struct of means of discrete approximations of
%                    null distributions; fieldnames match those of snr
%                    (listed above), except for 'sig'. Shapes and elements
%                    of arrays in each field also match and correspond one
%                    to one; e.g., the i_th, j_th element of
%                    snr.sig.null_mean.by_region is the mean of the null
%                    permutation distribution associated with the i_th,
%                    j_th element of snr.by_region.
%       .null_std  : Scalar struct of standard deviations of discrete
%                    approximations of null distributions; fieldnames match
%                    those of snr (listed above), except for 'sig'. Shapes
%                    and elements of arrays in each field also match and
%                    correspond one to one; e.g., the i_th, j_th element of
%                    snr.sig.null_mean.by_region is the standard deviation
%                    of the null permutation distribution associated with
%                    the i_th, j_th element of snr.by_region.
%
% Author: Jonathan Chien 5/22/22.


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addRequired(p, behavior);
addParameter(p, 'min_timeseries_length', 1); % Units of timepoints (25 Hz)
addParameter(p, 'bin_width', 10, @isscalar);
addParameter(p, 'step_size', 5, @isscalar);
addParameter(p, 'snr_dim', 1)
addParameter(p, 'zscore', true);
addParameter(p, 'normalization', 'bessel_corrected');
addParameter(p, 'decibels', true);
addParameter(p, 'pval', 'right-tailed');
addParameter(p, 'null_method', 'permutation')
addParameter(p, 'jitter_offset', 10);
addParameter(p, 'n_null', 1000);
addParameter(p, 'null_int_size', 95)
addParameter(p, 'return_null', false);
parse(p, timeseries, session_raw, behavior, varargin{:})
nv = p.Results;


%% Preprocessing

% Each region consists of a scalar timeseries.
n_regions = size(timeseries, 1);

% Pull out parts of timeseries labeled as behavior of interest.
timeseries_bhv = isolate_behavior(timeseries, session_raw, behavior);

% Pull out parts of timeseries labeled as baseline.
timeseries_base = isolate_baseline(timeseries, session_raw);

% Check that both baseline and behavior timeseries are of sufficient
% length. If timeseries are too short, they may for some regions consist
% only of zeros after deconvolution. Throw exception as caller if
% timeseries do not match or exceed min length.
if size(timeseries_bhv, 2) < nv.min_timeseries_length
   exception = MException("snr_by_bhv:timeseries_bhv_of_insufficient_length", ...
                          "The isolated timeseries for %s for current " + ...
                          "subject and session is shorter than the " + ...
                          "specified minimum length of %d (see " + ...
                          "'min_timeseries_length' name-value pair).", ...
                          behavior, nv.min_timeseries_length);
   throwAsCaller(exception);
elseif size(timeseries_base, 2) < nv.min_timeseries_length
    exception = MException("snr_by_bhv:timeseries_base_of_insufficient_length", ...
                           "The isolated baseline timeseries for current " + ...
                           "subject and session is shorter than the " + ...
                           "specified minimum length of %d (see " + ...
                           "'min_timeseries_length' name-value pair).", ...
                           nv.min_timeseries_length);
    throwAsCaller(exception);
end

% Bin data.
bin_timeseries_bhv = bin_data(timeseries_bhv, nv.bin_width, nv.step_size);
bin_timeseries_base = bin_data(timeseries_base, nv.bin_width, nv.step_size);


%% Calculate SNR by region

% This must be done before any z-scoring (else SNR will all be 1). Note
% that there will also not be any significance measures for this as
% permutation and jitter only affect covariance/correlation based metrics
% and will not change the computed SNR for a single region.
snr.by_region = snr_by_region(bin_timeseries_bhv, bin_timeseries_base, 0);


%% Calculate SNR upper bound

% Option to z-score data first. NB: if passed a vector of zeros, zscore
% will return all zeros. If a vector of zeros is passed to normalize.m, a
% vector of NaNs will be returned.
if nv.zscore
    % Z-score.
    bin_timeseries_bhv = zscore(bin_timeseries_bhv, 0, 2);
    bin_timeseries_base = zscore(bin_timeseries_base, 0, 2);
end

% Calculate overall SNR upper bound with all regions.
snr.full = filtered_snr(bin_timeseries_bhv, ...
                        bin_timeseries_base, ...
                        nv.snr_dim, ...
                        false, ... % Z-scoring already optionally performed above
                        nv.normalization, ...
                        nv.decibels);

% Calculate leave-one-out SNR upper bound.
snr.leave_one_out = leave_one_out_snr(bin_timeseries_bhv, ...
                                      bin_timeseries_base, ...
                                      n_regions, ...
                                      nv.snr_dim, ...
                                      nv.zscore, ...
                                      nv.normalization, ...
                                      nv.decibels);

% Compare leave-one-out SNR vs full SNR.
[snr.leave_one_out_diff, snr.leave_one_out_diff_normed] ...
    = leave_one_out_stats(snr.full, snr.leave_one_out);


%% Generate null distribution and attach p values.

if nv.pval
    % Generate null distributions. Unfortunately we need throwaway
    % variables that we can index into within the parfor loop because of
    % MATLAB's funny thing about not being able to prove indexing into
    % fields of structs within a parfor loop to be order independent.
    null_full = NaN(nv.snr_dim, nv.n_null);
    null_leave_one_out = NaN(n_regions, nv.snr_dim, nv.n_null);
    null_leave_one_out_diff = NaN(n_regions, nv.snr_dim, nv.n_null);
    null_leave_one_out_diff_normed = NaN(n_regions, nv.snr_dim, nv.n_null);

    parfor i_null = 1:nv.n_null
        % Generate one draw from underlying null.
        null_bin_timeseries_bhv ...
            = get_null_timeseries(bin_timeseries_bhv, ...
                                  nv.null_method, ...
                                  nv.jitter_offset);
        null_bin_timeseries_base ...
            = get_null_timeseries(bin_timeseries_base, ...
                                  nv.null_method, ...
                                  nv.jitter_offset);
    
        % Calculate SNR upper bound on null data.
        null_full(:,i_null) = filtered_snr(null_bin_timeseries_bhv, ...
                                           null_bin_timeseries_base, ...
                                           nv.snr_dim, ...
                                           nv.zscore, ...
                                           nv.normalization, ...
                                           nv.decibels);
    
        % Leave one out SNR upper bound on null data.
        null_leave_one_out(:,:,i_null) ...
            = leave_one_out_snr(null_bin_timeseries_bhv, ...
                                null_bin_timeseries_base, ...
                                n_regions, ...
                                nv.snr_dim, ...
                                nv.zscore, ...
                                nv.normalization, ...
                                nv.decibels);

        % Compare full vs leave-one-out.
        [null_leave_one_out_diff(:,:,i_null), ...
         null_leave_one_out_diff_normed(:,:,i_null)] ...
            = leave_one_out_stats(null_full(:,i_null), ...
                                  null_leave_one_out(:,:,i_null));
    end
    
    % Assign temporary vars into struct.
    null.full = null_full;
    null.leave_one_out = null_leave_one_out;
    null.leave_one_out_diff = null_leave_one_out_diff;
    null.leave_one_out_diff_normed = null_leave_one_out_diff_normed;
    
    % Compute significance related measures.
    fnames = fieldnames(snr);
    fnames(strcmp('by_region', fnames)) = []; % Remove 'by_region' fieldname as no sig metrics will be computed for it
    for i_field = 1:length(fnames)
        % Determine array dim containing null values.
        null_dim = length(size(null.(fnames{i_field}))); 

        % Attach p values.
        snr.sig.p.(fnames{i_field}) = tail_prob(snr.(fnames{i_field}), ...
                                                null.(fnames{i_field}), ...
                                                null_dim, ...
                                                'type', nv.pval, ...
                                                'exact', false, ...
                                                'correction', false);

        % 95 interval on null distribution.      
        snr.sig.null_int.(fnames{i_field}) ...
            = interval_bounds(null.(fnames{i_field}), ...
                              nv.null_int_size, ...
                              null_dim);

        % Calculate number of stdev from mean of null distribution, as well
        % as mean and std of null distribution.
        [snr.sig.z.(fnames{i_field}), ...
         snr.sig.null_mean.(fnames{i_field}), ...
         snr.sig.null_std.(fnames{i_field})] ...
            = null_zscore(snr.(fnames{i_field}), ...
                          null.(fnames{i_field}), ...
                          null_dim);
    end
end

end


% --------------------------------------------------
function snr = leave_one_out_snr(S, N, n_regions, n_dims, zscore_flag, ...
                                 normalization, decibel_flag)
% Calculate leave-one-out SNR excluding in turn each of the regions. Output
% is an m x n array, where m = number of regions and n = dimensionality of
% SNR.

snr = NaN(n_regions, n_dims); 
for i_region = 1:n_regions
    snr(i_region,:) = filtered_snr(S(setdiff(1:n_regions, i_region),:), ...
                                   N(setdiff(1:n_regions, i_region),:), ...
                                   n_dims, zscore_flag, normalization, ...
                                   decibel_flag);
end

end


% --------------------------------------------------
function [delta, delta_normed] = leave_one_out_stats(full, leave_one_out)
% Compare SNR with each region left out vs SNR with all regions. Typically,
% full is of size n x 1, and leave_one_out is of size m x n, where n = the
% number of dimensions of the SNR and m = the number of regions.

delta = full' - leave_one_out;

delta_normed = delta ./ full';

end

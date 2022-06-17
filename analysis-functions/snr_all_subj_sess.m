function [snr_out_by_bhv, snr_out_by_subj_sess, rmv_bhv_ind] ...
    = snr_all_subj_sess(path_to_data, behavior_list, subj_ind, varargin)
% Accepts as input a list of behaviors and subject index range. For each
% subject, this function attempts to calculate the SNR upper bound, leave
% one out SNR upper bound, and by-region SNR for each behavior in each
% session (not all behaviors occur in all sessions). All significance
% measures are computed by first either randomly permuting or jittering the
% timeseries and then repeating all calculations, for each independent
% random permutation, on the permuted data.
% 
% PARAMETERS
% ----------
% path_to_data  : Path to directory containing data.
% behavior_list : n_behaviors x 1 cell array whose i_th cell contains the
%                 name of the i_th behavior to be attempted.
% subj_ind      : n_subjects x 1 array of indices corresponding to
%                 subjects. If one of the indices does not correspond to a
%                 filename, the function will skip it. 
% Name-Value Pairs (nv)
%   'data_type'          : ('raw' (default) | 'devonv' | 'mpp'). Set to
%                          'raw' to use data as passed in. Set to 'deconv'
%                          to apply deconvolution via the OASIS package
%                          (https://github.com/zhoupc/OASIS_matlab). Set to
%                          'mpp' to use peak-finding/convolution for
%                          conversion of data to a filtered Markov point
%                          process. Note that visualization of the
%                          deconvolved traces for one subject/session is
%                          supported in the local function wrapping
%                          deconvolveCa.m, apply_deconvolution.
%   'preprocess_params'  : Scalar struct whose fields contain the
%                          parameters for preprocessing of data (ignored if
%                          'data_type' = 'raw'). The fields should
%                          correspond to optional arguments of the function
%                          deconvolveCa.m for 'deconv'
%                          markov_point_process for 'mpp'.
%   'snr_params'         : Scalar struct whose fields contain the
%                          name-value pair arguments for snr_by_bhv.m
%   'min_session_number' : (Scalar), specify the minimum number of sessions
%                          in which a behavior needs to appear for that
%                          behavior to be retained.
%
% RETURNS
% snr_out_by_bhv       : Scalar struct. Each field is named for a behavior
%                        and contains an n_sessions_i x 1 struct, where
%                        n_sessions_i is the number of sessions, across
%                        subjects, in which the given behavior occurred.
%                        The i_th element of each of the fields of the
%                        n_sessions_i x 1 struct are, for the i_th session
%                        featuring the given behavior:
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
%                                field of snr_params is not false and has a
%                                valid value. If valid, sig is a scalar
%                                struct with the following fields:
%       .p         : Scalar struct of p values whose fieldnames match those
%                    listed above (for one session from each behavior)
%                    under snr_out_by_bhv (except for sig). Array shapes
%                    also match and elements correspond one to one; e.g.,
%                    the i_th, j_th element of .sig.p.by_region is the p
%                    value associated with the i_th, j_th element of
%                    .by_region.
%       .null_int  : Scalar struct of null intervals whose fieldnames match
%                    those listed above (for one session from each
%                    behavior) under snr_out_by_bhv (except for sig). For
%                    each field under one behavior/session in
%                    snr_out_by_bhv, the corresponding field in
%                    sig.null_int has an array with one extra dimension (to
%                    store the lower and upper interval bounds). For
%                    example, if .leave_one_out has shape 13 x 2, then
%                    sig.null_int_leave_one_out has shape 13 x 2 x 2, where
%                    the (i,j,:) slice contains the upper and lower
%                    interval bounds for the i_th, j_th value of
%                    .leave_one_out.
%       .z         : Scalar struct of z scores whose fieldnames match those
%                    listed above (for one session from each behavior)
%                    under snr_out_by_bhv (except for sig). Array shapes
%                    also match and elements correspond one to one; e.g.,
%                    the i_th, j_th element of .sig.z.by_region is the
%                    z-score on the null permutation distribution
%                    associated with the i_th, j_th element of .by_region.
%       .null_mean : Scalar struct of means of discrete approximations of
%                    null distributions; fieldnames match those listed
%                    above (for one session from each behavior) under
%                    snr_out_by_bhv (except for sig). Array shapes also
%                    match and elements correspond one to one; e.g., the
%                    i_th, j_th element of .sig.null_mean.by_region is the
%                    mean of the null permutation distribution associated
%                    with the i_th, j_th element of .by_region.
%       .null_std  : Scalar struct of standard deviations of discrete
%                    approximations of null distributions; fieldnames match
%                    those listed above (for one session from each
%                    behavior) under snr_out_by_bhv (except for sig). Array
%                    shapes also match and elements correspond one to one;
%                    e.g., the i_th, j_th element of
%                    .sig.null_std.by_region is the standard deviation of
%                    the null permutation distribution associated with the
%                    i_th, j_th element of .by_region.
% snr_out_by_subj_sess : n_subj_ind x n_sessions_max struct, where
%                        n_subj_ind is the number of elements passed in
%                        through the subj_ind argument (see PARAMETERS) and
%                        n_sessions_max is the most number of sessions any
%                        subject had. The i_th, j_th element contains, for
%                        the j_th session of the i_th subject, a
%                        n_behaviors x 1 struct, where each field contains
%                        the same values as one element of the fields of
%                        the nonscalar structs corresponding to one
%                        behavior and session in snr_out_by_bhv (these
%                        fields are by_region, full, leave_one_out,
%                        leave_one_out_diff, leave_one_out_diff_normed).
% rmv_bhv_ind          : b x 1 vector of behavior indices, where b is the
%                        number of behaviors omitted. After computation of
%                        SNRs, the user will be warned if any behaviors
%                        were skipped due to insufficient timeseries length
%                        (see 'min_timeseries_length' name-value pair in
%                        snr_by_bhv) or due to an insufficient number of
%                        sessions (see 'min_session_number' above). The
%                        indices of any behaviors (according to their
%                        position in the behavior_list argument (see
%                        PARAMETERS)) will be returned here.
%
% Author: Jonathan Chien. 5/24/22.


%% Parse inputs

p = inputParser;
addRequired(p, 'path_to_data');
addRequired(p, 'behavior_list');
addRequired(p, 'subj_ind');
addParameter(p, 'data_type', 'raw');
addParameter(p, 'preprocess_params', [])
addParameter(p, 'snr_params', []);
addParameter(p, 'min_session_number', 1) % Minimum number of sessions each behavior must have, else it is excluded
parse(p, path_to_data, behavior_list, subj_ind, varargin{:})
nv = p.Results;


%% Setup

% Get array dim sizes and preallocate.
n_subjects_to_try = length(subj_ind);
n_behaviors_to_try = length(behavior_list);
MAX_N_SESSIONS = 100; % A number which no subject's session count will exceed
snr = cell(n_subjects_to_try, MAX_N_SESSIONS, n_behaviors_to_try);


%% Calculate SNR upper bound for all behaviors over all subjects/sessions

for i_subj = subj_ind
    % Load all sessions for single subject.
    try
        [sessions_raw, n_regions_retained, nan_regions, session_ids] ...
            = extract_single_subj(path_to_data, i_subj, ...
                                  'selection_ind', 'all', ...
                                  'remove_nan', true, ...
                                  'nan_warning', false, ...
                                  'plot', false);
        n_sessions = length(sessions_raw);
    catch exception
        % Allow continuation of program only if error is as expected due to
        % no subject having the current index. Else, rethrow/terminate.
        if strcmp(exception.identifier, "extract_single_subject:file_not_found")
            continue 
        else
            rethrow(exception);
        end
    end

    % Initialize waitbar.
    w = waitbar(0, '');

    % Test all behaviors possible in all sessions for current subject.
    for i_sess = 1:n_sessions
        %---------------------Optional preprocessing ---------------------%

        switch nv.data_type
            % Deconvolution.
            case 'deconv'
                assert(~isempty(nv.preprocess_params), ...
                       "If 'preprocess' is not false, you must pass in the " + ...
                       "requisite parameters in the fields of a scalar " + ...
                       "struct provided as the value for the preprocess_params " + ...
                       "name-value pair.")
                waitbar((i_sess-1)/n_sessions, w, ...
                         sprintf("Deconvolving data for subject %d: " + ...
                                 "session %d of %d...", ...
                                 i_subj, i_sess, n_sessions));
                neural_data ...
                    = apply_deconvolution(sessions_raw{i_sess}, ...
                                          nv.preprocess_params, ...
                                          false); % Plotting
            % Conversion to filtered Markov point process.
            case 'mpp'
                assert(~isempty(nv.preprocess_params), ...
                       "If 'preprocess' is not false, you must pass in the " + ...
                       "requisite parameters in the fields of a scalar " + ...
                       "struct provided as the value for the preprocess_params " + ...
                       "name-value pair.")
                waitbar((i_sess-1)/n_sessions, w, ...
                        sprintf("Converting subject %d: session %d of %d " + ...
                                "data to filtered Markov point process...", ...
                                i_subj, i_sess, n_sessions))
                neural_data ...
                    = markov_point_process(sessions_raw{i_sess}.Lfold, ...
                                           'min_threshold', nv.preprocess_params.min_threshold, ...
                                           'zscore', nv.preprocess_params.zscore, ...
                                           'convolve', nv.preprocess_params.convolve, ...
                                           'plot', false);
            case 'raw'
                neural_data = sessions_raw{i_sess}.Lfold;
            otherwise
                error("Please provide a valid value for 'data_type'.")
        end

        % Update waitbar.
        waitbar(i_sess/n_sessions, w, ...
                sprintf("Working on subject %d: session %d of %d...", ...
                        i_subj, i_sess, n_sessions));


        % -------------------------Calculate SNR -------------------------%
        
        % Ensure that the requested number of SNR dimensions does not
        % exceed the number of regions retained. Throw error if so.
        if nv.snr_params.snr_dim > n_regions_retained(i_sess)
            error("The number of requested SNR dimensions exceeds the " + ...
                  "number of regions retained for subject %d session %s. " + ...
                  "Set snr_params.snr_dim to a lower value.", ...
                  i_subj, session_ids{i_sess});
        end

        for i_bhv = 1:n_behaviors_to_try
            % Ensure user passed in name-value pair params for snr_by_bhv.
            assert(~isempty(nv.snr_params), ...
                   "You must pass in name-value pairs for snr_wrapper.m " + ...
                   "as the fields of the struct passed as value for the " + ...
                   "name 'snr_params'.")

            % Calculate optimized SNR.
            try
                snr{i_subj,i_sess,i_bhv} = snr_by_bhv( ...
                        neural_data, ...
                        sessions_raw{i_sess}, ...
                        behavior_list{i_bhv}, ...
                        'min_timeseries_length', nv.snr_params.min_timeseries_length, ...
                        'bin_width', nv.snr_params.bin_width, ...
                        'step_size', nv.snr_params.step_size, ...
                        'zscore', nv.snr_params.zscore, ...
                        'normalization', nv.snr_params.normalization, ...
                        'decibels', nv.snr_params.decibels, ...
                        'snr_dim', nv.snr_params.snr_dim, ...
                        'pval', nv.snr_params.pval, ...
                        'n_null', nv.snr_params.n_null, ...
                        'null_int_size', nv.snr_params.null_int_size, ...
                        'return_null', nv.snr_params.return_null ...
                                                      );
            catch exception
                if strcmp(exception.identifier, ...
                          "isolate_behavior:behavior_not_found")
                    % Allow continuation on to next behavior if exception
                    % is due to behavior not being found in current
                    % session.
                    continue
                elseif strcmp(exception.identifier, ...
                              "isolate_baseline:unable_to_isolate_baseline")
                    % Skip entire current session if a baseline could not
                    % be identified. Issue warning if this happens.
                    warning("Session %s for subject %d will be skipped " + ...
                            "(for all behaviors) since a baseline period " + ...
                            "could not be identified", ...
                            session_ids{i_sess}, i_subj)
                    break
                elseif strcmp(exception.identifier, ...
                              "snr_by_bhv:timeseries_bhv_of_insufficient_length") ...
                       || strcmp(exception.identifier, ...
                                 "snr_by_bhv:timeseries_base_of_insufficient_length")
                    % Allow continuation on to next behavior if exception
                    % is due to insufficient length of extracted
                    % timeseries.
                    continue          
                else
                    % If exception is not due to the above, rethrow and
                    % terminate for debugging.
                    rethrow(exception)
                end
            end

            
            % ----------------------NaN reinsertion ----------------------%

            % If regions were removed due to NaN values, add the NaN
            % values back here so that the region index consistently
            % matches the same region across all subjects/sessions. 
            if ~isempty(nan_regions{i_sess})
                % First determine total number of regions for current
                % subject/session.
                n_regions_tot = n_regions_retained(i_sess) ...
                                + length(nan_regions{i_sess});
                
                % Get names of fields containing arrays that may require
                % NaN reinsertion.
                may_need_mod = fieldnames(snr{i_subj,i_sess,i_bhv});
                for i_field = 1:length(may_need_mod)
                    % If array is not from the 'full' field, we will need
                    % to re-insert NaN values.
                    if ~ismember(may_need_mod{i_field}, {'full', 'sig'})
                        snr{i_subj,i_sess,i_bhv}.(may_need_mod{i_field}) ...
                            = reinsert_nan(snr{i_subj,i_sess,i_bhv}.(may_need_mod{i_field}), ...
                                           n_regions_tot, ...
                                           nan_regions{i_sess});

                        % Also reinsert NaN value into sigifiance measures 
                        % if they were computed.
                        if isfield(snr{i_subj,i_sess,i_bhv}, 'sig') ...
                           && ~strcmp(may_need_mod{i_field}, 'by_region') % NaNs need to be reinserted for by_region SNR, but it has no corresponding sig metrics
                            sig_measures = fieldnames(snr{i_subj,i_sess,i_bhv}.sig);
                            for ii_field = 1:length(sig_measures)                               
                                snr{i_subj,i_sess,i_bhv} ...
                                    .sig ...
                                    .(sig_measures{ii_field}) ...
                                    .(may_need_mod{i_field}) ...
                                    = reinsert_nan(snr{i_subj,i_sess,i_bhv} ...
                                                       .sig ...
                                                       .(sig_measures{ii_field}) ...
                                                       .(may_need_mod{i_field}), ...
                                                   n_regions_tot, ...
                                                   nan_regions{i_sess});                           
                            end
                        end
                    end        
                end
            end
        end    
    end

    close(w);
end


%% Clean up

% First mark for removal any behaviors for which we did not/were unable to
% calculate SNR (for example, if all sessions for that behavior did not
% match or exceed the specified minimum timeseries length--see
% 'min_timeseries_length' name-value pair for snr_by_bhv.m). Issue warning
% to user if any behaviors were marked for removal.
empty_bhv_ind = squeeze(all(cellfun(@isempty, snr), [1 2]));
if any(empty_bhv_ind)
    warning("All sessions for %s will be removed, most likely due to no " + ...
            "sessions for this/these behavior(s) matching or exceeding " + ...
            "the specified minimum timeseries length (see " + ...
            "'min_timseries_length' name-value-pair in snr_by_bhv).", ...
            strjoin(behavior_list(empty_bhv_ind), ", "))
end

% Next, mark for removal any behaviors that did not have enough sessions
% (see 'min_session_number' name-value pair from this function). Issue
% warning if any behaviors identified for removal.
insuff_session_ind = squeeze(sum(~cellfun(@isempty, snr), [1 2])) ...
                     < nv.min_session_number;
if any(insuff_session_ind)
    warning("Session(s) for %s will be removed, as the number of sessions " + ...
            "was less than the minimum (%d) specified through the " + ...
            "'min_session_number' name-value pair.", ...
            strjoin(behavior_list(insuff_session_ind), ", "), ...
            nv.min_session_number)
end

% Get indices of all behaviors to be removed and remove those behaviors
% from both container array and list of behaviors.
rmv_bhv_ind = empty_bhv_ind | insuff_session_ind;
snr(:,:,rmv_bhv_ind) = [];
behavior_list(rmv_bhv_ind) = [];

% Next remove extra cells from 2nd dim of cell array (consisting of extra
% padding for number of sessions).
max_n_sessions = max(sum(~cellfun(@isempty, snr), 2), [], 'all');
snr = snr(:,1:max_n_sessions,:);

% Associate each slice along third cell array dim with corresponding
% behavior by creating scalar struct with fields corresponding to behaviors
% and each such field containing an n_subject x max_n_sessions cell array
% where each cell contains the scalar struct output of snr_by_bhv.
snr_out_by_bhv = struct();
for i_bhv = 1:length(behavior_list)
    snr_out_by_bhv.(behavior_list{i_bhv}) = squeeze(snr(:,:,i_bhv));
end

% Create an n_subjects x n_behaviors cell array where each cell contains a
% struct with fields corresponding to behaviors, with each such field being
% the scalar struct output by snr_by_bhv.
snr_out_by_subj_sess = cell2struct(snr, behavior_list, 3);

end


% --------------------------------------------------
function deconv = apply_deconvolution(session_raw, params, plot_flag)
% Wrapper for deconvolution. params corresponds to nv.preprocess_params in
% the invoking host function.

% Preallocate.
deconv = NaN(size(session_raw.Lfold));
n_regions = size(deconv, 1);
        
for i_region = 1:n_regions
    % Deconvolve.
    [~, deconv(i_region,:), ~] ...
        = deconvolveCa(session_raw.Lfold(i_region,:), ...
                       params.kernel_model, params.kernel_params, ...
                       params.method, params.extra_params, ...
                       'window', params.window, ...
                       'sn', params.noise_std, ...
                       'b', params.baseline, ...
                       'lambda', params.lambda, ...
                       'optimize_b', params.optimize_b, ...
                       'optimize_pars', params.optimize_pars, ...
                       'optimize_smin', params.optimize_smin);       
end

% Optional plotting. Helpful for debugging/pausing for visualization.
if plot_flag
    for i_region = 1:n_regions
        % Prepare window.
        if true
            plot_window = [1 size(deconv, 2)];
        else
            plot_window = [2000 5000]; % Set desired window
        end
        
        % Plot deconvolved.
        subplot(n_regions, 1, i_region)
        hold on
        stem(deconv(i_region,:), 'Marker', 'none')
        xlim(plot_window)
    end
end

end


% --------------------------------------------------
function Y = reinsert_nan(X, n_regions_tot, nan_regions)
% Currently, n (where n = number of regions after NaN removal) is the 1st
% array dim size. This convention must not be changed. An assertion clause
% adds some added protection but this construction is admittedly not the
% most robust. X is a nonscalar array whose 1st dim size is n and whose
% other dims may correspond to dimensionality of SNR signal, upper or lower
% interval bounds etc. Returned is Y, an array whose 1st dim size is m,
% where m = n + n_regions_removed, with all other dim sizes matching those
% of X; the added entries consist of NaN padding.

assert(size(X, 1) + length(nan_regions) == n_regions_tot);

sizes = size(X); sizes(1) = n_regions_tot;
Y = zeros(sizes);
Y(nan_regions,:,:) = NaN; % Recall that extra trailing colons are allowed/ignored
Y(~isnan(Y)) = X;

end

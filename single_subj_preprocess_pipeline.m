
clearvars
close all
clc


%% Load data for one subject

% Specify data location and subject to be analyzed.
if strcmp(computer, 'PCWIN64')
    path_to_data = '';
else
    path_to_data = '';
end
i_subject = 1;
selected_session_ind = 'all'; % 'all' or vector of indices

% Set plotting parameters.
params.plotting.raw.plot = false;
params.plotting.raw.figure = true;
params.plotting.raw.window = [1000 3000];
params.plotting.raw.scale_factor = 0.3;

% Extract all sessions (or subset of sessions) data for one subject.
[sessions_raw, n_regions, nan_regions] ...
    = extract_single_subj(path_to_data, i_subject, ...
                          'selection_ind', [1 2], ...
                          'plot', params.plotting.raw.plot, ...
                          'figure', params.plotting.raw.figure, ...
                          'window', params.plotting.raw.window, ...
                          'scale_factor', params.plotting.raw.scale_factor);
                      
% Get number of sessions extracted.
n_selected_sessions = length(sessions_raw);


%% Cross-correlation analysis for regions by behavior (raw traces)

% close all 
clc

% Specify session (among those extracted) to be analyzed here.
i_session = 2;
behavior = 'Investigate';

params.xcorr.raw.bin_width = 20;
params.xcorr.raw.step_size = 5;
params.xcorr.raw.zscore = true;
params.xcorr.raw.scaleopt = 'unbiased';
params.xcorr.raw.maxlag = 10;
params.xcorr.raw.svn_opts = [false false]; % Options to normalize by baseline and scale down by std after calculating snr
params.xcorr.raw.pval = false;
params.xcorr.raw.n_null = 1000;
params.xcorr.raw.null_method = 'jitter';
params.xcorr.raw.jitter_offset = 20;
params.xcorr.raw.threshold = 0.05;
params.xcorr.raw.extract_opt_by = 'max_xcorr';
params.xcorr.raw.cmap = parula;

% Exploratory. This function will allow visualization of cross-correlations
% among all regions across a range of lags, with p values attaching to all
% cross-correlations via MC permutations.
x_corr_raw = xcorr_by_bhv(sessions_raw{i_session}.Lfold, ...
                          sessions_raw{i_session}, ...
                          behavior, ...
                          'bin_width', params.xcorr.raw.bin_width, ...
                          'step_size', params.xcorr.raw.step_size, ...
                          'zscore', params.xcorr.raw.zscore, ...
                          'scaleopt', params.xcorr.raw.scaleopt, ...
                          'maxlag', params.xcorr.raw.maxlag, ...
                          'svn_opts', params.xcorr.raw.svn_opts, ...
                          'pval', params.xcorr.raw.pval, ...
                          'n_null', params.xcorr.raw.n_null, ...
                          'null_method', params.xcorr.raw.null_method, ...
                          'return_null', true, ...
                          'threshold', params.xcorr.raw.threshold, ...
                          'extract_opt_by', params.xcorr.raw.extract_opt_by, ...
                          'plot_all_lags', true, ...
                          'plot_opt_lag', 'corr_mat', ...
                          'cmap', params.xcorr.raw.cmap);
                               
 
 %% Optional deconvolution

% Specify whether or not to deconvolve and if so, whether or not to
% visualize results (if plot desired, set window width either as 2-vector
% or as the string 'entire').
params.deconv.convolve = true;

% Set parameters related to kernel.
params.deconv.kernel_model = 'ar1'; 
params.deconv.kernel_params = []; % Corresponds to 'pars', default = []
params.deconv.method = 'constrained'; 
params.deconv.extra_params = []; % Only if method = 'mcmc'
params.deconv.window = 200; % Width of kernel

% Other parameters.
params.deconv.noise_std = []; % Default is [], in which case this will be estimated by function
params.deconv.baseline = 0; % Default = 0
params.deconv.lambda = 0; % Default = 0
params.deconv.optimize_b = false;
params.deconv.optimize_pars = false;
params.deconv.optimize_smin = false;

% As of 5/19/22, any parameters not specified should be assumed to match
% the defaults at the time the repo was cloned on 5/18/22.

% Preallocate for sessions.
deconv = cell(n_selected_sessions, 1);

for i_session = 1:n_selected_sessions
    % Preallocate witihin session.
    deconv{i_session} = NaN(size(sessions_raw{i_session}.Lfold));
        
    for i_region = 1:n_regions(i_session) 
        % Deconvolve.
        [~,deconv{i_session}(i_region,:),~] ...
            = deconvolveCa(sessions_raw{i_session}.Lfold(i_region,:), ...
                           params.deconv.kernel_model, params.deconv.kernel_params, ...
                           params.deconv.method, params.deconv.extra_params, ...
                           'window', params.deconv.window, ...
                           'sn', params.deconv.noise_std, ...
                           'b', params.deconv.baseline, ...
                           'lambda', params.deconv.lambda, ...
                           'optimize_b', params.deconv.optimize_b, ...
                           'optimize_pars', params.deconv.optimize_pars, ...
                           'optimize_smin', params.deconv.optimize_smin);       
    end
end 


%% Option to visualize deconvolved data (if deconvolution performed).

params.plotting.deconv.plot = true;
params.plotting.deconv.window = [500 5000];

if params.conv.convolve && params.plotting.deconv.plot
    figure
    i_subplot = 1;
    
    for i_session = 1:n_selected_sessions
    for i_region = 1:n_regions(i_session)
        % Prepare window.
        if strcmp(params.plotting.deconv.window, 'entire')
            plot_window = [1 length(deconv{i_session})];
        else
            plot_window = params.plotting.deconv.window;
        end
        
        % Plot current session.
        subplot(n_regions(i_session), n_selected_sessions, i_subplot)
        hold on
        stem(deconv{i_session}(i_region,:), 'Marker', 'none')
        xlim(plot_window)
%         xlabel(sprintf('Time at 25 Hz from %dth to %dth timepoint', plot_window(1), plot_window(2)))
%         title(sprintf('Session %s', selected_session_ids{i_session}))
        clear plot_window
        
        i_subplot = i_subplot + 1;
    end  
    end
end


%% Cross-correlations on deconvolved data

% Specify session and behavior to be analyzed.
i_session = 1;
behavior = 'Investigate';

params.xcorr.deconv.bin_width = 6;
params.xcorr.deconv.step_size = 2;
params.xcorr.deconv.zscore = true;
params.xcorr.deconv.scaleopt = 'unbiased';
params.xcorr.deconv.maxlag = 10;
params.xcorr.raw.snr_opts = [true true];
params.xcorr.deconv.pval = 'two-tailed';
params.xcorr.deconv.n_null = 800;
params.xcorr.deconv.null_method = 'jitter';
params.xcorr.deconv.threshold = false;

x_corr_deconv = xcorr_by_bhv(deconv{i_session}, ...
                             sessions_raw{i_session}, ...
                             behavior, ...
                             'bin_width', params.xcorr.deconv.bin_width, ...
                             'step_size', params.xcorr.deconv.step_size, ...
                             'zscore', params.xcorr.deconv.zscore, ...
                             'scaleopt', params.xcorr.deconv.scaleopt, ...
                             'maxlag', params.xcorr.deconv.maxlag, ...
                             'snr_opts', params.xcorr.raw.snr_opts, ...
                             'pval', params.xcorr.deconv.pval, ...
                             'n_null', params.xcorr.deconv.n_null, ...
                             'null_method', params.xcorr.deconv.null_method, ...
                             'return_null', false, ...
                             'threshold', params.xcorr.deconv.threshold, ...
                             'plot_flattened', true, ...
                             'plot_max', true, ...
                             'cmap', 'parula');


%% Option to use peak detection + convolution to create filtered Markov point process

% Set threshold (as percentage of max activity for each region) below which
% all activity will be set to zero).
threshold_percent = 0.25;
filter = [0.1 0.2 0.4 0 0] / 0.7;

figure

mpp = cell(n_selected_sessions, 1);
for i_session = 1:n_selected_sessions
    mpp{i_session} = markov_point_process(sessions_raw{i_session}.Lfold, ...
                                          'min_threshold', threshold_percent, ...
                                          'normalize', true, ...
                                          'convolve', filter, ...
                                          'plot', true, ...
                                          'figure', false, ...
                                          'window', [3000 7500]);
end


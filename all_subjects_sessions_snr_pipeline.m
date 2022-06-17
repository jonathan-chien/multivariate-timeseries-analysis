
clearvars
clc
close all


%% Get list of all behaviors

path_to_data = '';
all_behaviors = get_all_behaviors(path_to_data, [1 33]);


%% Load all subjects/sessions and calculate SNR upper bound for all behaviors

% Set whether or not to flatten outputs (necessary for subsequent plotting
% steps).
FLATTEN = true;

% Specify which behaviors to exclude.
excluded_behaviors = {'Intro_Baseline', 'Rmv_Baseline', 'Other'};
behaviors_to_test = setdiff(all_behaviors, excluded_behaviors);
subj_ind = 1:33; % Indices are for subjects


%------------------------ Preprocessing parameters -----------------------%

% Specify what data type to use.
params.data_type = 'raw'; % {'raw' | 'deconv' | 'mpp'}; note that 'mpp' is not recommended

% This is just a placeholder to ensure there is a valid field name during
% call to snr_by_bhv, as there are no hyperparameters for preprocessing if
% we are using raw traces.
params.raw = [];

% Set parameters for possible deconvolution.
params.deconv.kernel_model = 'ar1'; 
params.deconv.kernel_params = []; % Corresponds to 'pars', default = []
params.deconv.method = 'constrained'; 
params.deconv.extra_params = []; % Only if method = 'mcmc'
params.deconv.window = 8; % Width of kernel
params.deconv.noise_std = []; % Default is [], in which case this will be estimated by function
params.deconv.baseline = 0; % Default = 0
params.deconv.lambda = 0; % Default = 0
params.deconv.optimize_b = false;
params.deconv.optimize_pars = false;
params.deconv.optimize_smin = false;

% Set parameters for possible peak detection/filtering.
params.mpp.min_threshold = 0.25;
params.mpp.zscore = true;
params.mpp.convolve = [0.1 0.2 0.4 0 0] / 0.7;


%----------------------------- SNR parameters ----------------------------%

params.snr.min_timeseries_length = 50; % Units of timepoints (Hz); 100 = 4 s
params.snr.min_session_number = 5; % Minimum number of sessions for a behavior to be retained; note this is passed directly to snr_all_subjects_sessions 
params.snr.bin_width = 8;
params.snr.step_size = 4;
params.snr.zscore = true;
params.snr.normalization = 'trace';
params.snr.decibels = true;
params.snr.snr_dim = 2;
params.snr.pval = 'right-tailed';
params.snr.n_null = 8;
params.snr.null_int_size = 95;
params.snr.return_null = true;


%------------- Calculate SNR upper bound for each behavior----------------%

% SNR for all behaviors in all sessions.
[snr_out_by_bhv, snr_out_by_subj_sess, rmv_bhv_ind] ...
    = snr_all_subj_sess(path_to_data, behaviors_to_test, subj_ind, ...
                        'data_type', params.data_type, ...
                        'preprocess_params', params.(params.data_type), ...
                        'snr_params', params.snr, ...
                        'min_session_number', params.snr.min_session_number);
behaviors_to_test(rmv_bhv_ind) = [];

% For each behavior, flatten cell array into struct array.
if FLATTEN
    for i_bhv = 1:length(behaviors_to_test)  
        snr_out_by_bhv.(behaviors_to_test{i_bhv}) ...
            = vertcat(snr_out_by_bhv.(behaviors_to_test{i_bhv}){:});
    end
end


%% Plot data 

clc
close all

% Plot each session's SNR 1 and 2 for each behavior.
figure
plot_session_snr(snr_out_by_bhv, behaviors_to_test, "sig_thresh", 0.05)

% Plot one grid for each behavior; each grid is n_subjects x n_sessions and
% displays SNR 1.
figure
plot_snr_by_bhv_subj_sess(snr_out_by_subj_sess, behaviors_to_test);
sgtitle('Z-scored SNR')

% Plot each behavior's mean SNR 1 across sessions.
figure
plot_mean_snr(snr_out_by_bhv, behaviors_to_test, 1, 'convert_to_decibel', true)



%% PCA/SVD analysis of behavior-region interactions

close all
clc

% Load region names.
load("")

% Extract data and perform PCA.
[pca_out, pca_array] = snr_pca(snr_out_by_bhv, ...
                               "observations", "single_session", ...
                               "normalize", true);

% Plot biplot.
figure
biplot(pca_out.loadings(:,1:3), 'VarLabels', regions)

% Plot eigenspectrum.
figure
plot(pca_out.eigenspectrum, 'LineWidth', 3);
xlabel('PC')
ylabel('Eigenvalue')
title('Eigenspectrum')
xlim([1 length(pca_out.eigenspectrum)])


% First set number of PCs to plot. 
n_pcs_to_plot = 3; pc_labels = cell(n_pcs_to_plot, 1);

% Generate labels for PCs.
for k = 1:n_pcs_to_plot, pc_labels{k} = sprintf("PC %d", k); end

% Generate parallel coordinate plot.
figure
parallel_coord_plot(pca_out.U(:,1:n_pcs_to_plot)', pc_labels, [], ...
                     'same_plot', false)


function [sessions_raw, n_regions_retained, nan_regions, selected_session_ids] ...
    = extract_single_subj(path_name, i_subject, varargin)
% Extracts all raw data from all sessions (or subset of sessions) for one
% subject. Optional visualization of raw traces.
%
% PARAMETERS
% ----------
% path      : String, path to the directory containing raw data.
% i_subject : Scalar integer, index of the subject to be analyzed.
% Name-Value Pairs
%   'selection_ind' : ('all' (default) | vector), set as 'all' to load
%                      data from all sessions. Else, pass in a vector of
%                      indices of sessions we wish to load; indices must be
%                      from the set {1:n}, where n = the number of sessions
%                      for the current subject. The length of this vector
%                      is equal to n_sessions_selected (see RETURNS).
%   'remove_nan'    : (1 (default) | 0), specify whether to remove regions
%                     that have NaN values.
%   'nan_warning'   : (1 (default) | 0), specify whether to issue warning
%                     if regions were removed due to NaN values. Ignored if 
%                     'remove_nan' = false.
%   'plot'          : (1 (default) | 0), specify whether or not to
%                      visualize raw traces.
%   'figure'        : (1 (default) | 0), specify whether or not to plot
%                     into new figure. Ignored if 'plot' is false.
%   'scale_factor'  : Scalar value multiplying the raw traces; can be used
%                     to scale down traces to prevent overlap. Ignored if
%                     'plot' is false.
%   'window'        : 2-vector whose 1st and 2nd elements are the indices
%                     of the 1st and 2nd timepoints, respectively, of the
%                     window we would like to visualize. Ignored if 'plot'
%                     is false.
%
% RETURNS
% -------
% sessions_raw         : n_sessions_selected x 1 cell array. The k_th cell
%                        contains a scalar struct with the following
%                        fields:
%   .Lfold     : Recorded GCaMP traces
%   .t         : Time stamp of each GCaMP imaging frame
%   .behaviors : Annotated behavior
%   .Fstart    : Starting frame of behaviors
%   .FStop     : Stopping frame of behaviors
%   .regions   : Brain region
%   .FL        : Sampling rate, 25 Hz --> 40 ms
% n_regions_retained   : n_sessions_selected x 1 array whose i_th element
%                        is the number of regions retained after possible
%                        NaN removal for the i_th session.
% nan_regions          : This is an n_sessions_selected x 1 cell array
%                        whose i_th cell contains a vector whose elements
%                        are the indices of regions removed due to
%                        containing NaN values. If 'nan_regions' = false,
%                        all cells will be empty.
% selected_session_ids : n_sessions_selected x 1 cell array whose i_th cell
%                        contains the session ID of the i_th session
%                        selected/to be returned.
%
% Author: Jonathan Chien. 5/19/22. Refactored from
% preprocess_single_subject_multi_sessions.


%% Parse inputs

p = inputParser;
addRequired(p, 'path', @ischar);
addRequired(p, 'i_subject', @(x) isscalar(x) && all(floor(x) == x));
addParameter(p, 'selection_ind', 'all', @(x) ischar(x) || all(floor(x) == x))
addParameter(p, 'remove_nan', true, @islogical)
addParameter(p, 'nan_warning', true, @islogical)
addParameter(p, 'plot', true, @(x) islogical(x))
addParameter(p, 'figure', true, @(x) islogical(x))
addParameter(p, 'scale_factor', 0.3, @(x) isscalar(x))
addParameter(p, 'window', [1000 3000], @(x) length(x) == 2 && all(floor(x) == x))
parse(p, path_name, i_subject, varargin{:});
nv = p.Results;


%% Load data

% Add directory containing data to path.
c = pathsep;
path_str = [c path c];
on_path = contains(path_str, [c path_name c], 'IgnoreCase', ispc);
if ~on_path, addpath(path_name); end

% Set subject index and load data. Since there are technically many ways to
% generate an error with load, this setup here probably isn't a great way
% to do things, but we will treat any errors as due to the filename not
% existing and throw an exception in invoking function named (as of 6/1/22)
% snr_all_subjects_session; this will at least cause errors within this
% function not related to file-loading to be caught.
filename = sprintf('MFP-ERa-GC%d.mat', i_subject);
try
    data_as_loaded = load(filename);
catch
    exception = MException("extract_single_subject:file_not_found", ...
                           "Filename %s not found", filename);
    throwAsCaller(exception);
end

% Determine number of sessions and get session IDs.
fnames = fieldnames(data_as_loaded.Raw);
n_sessions = length(fnames) - 6; % There are six fields before the sessions
all_session_ids = fnames(7 : end);


%% Select data from all sessions or from subset of sessions

% Set the indices of the sessions we would like to run. Set 'all' or empty
% array to load all sessions.
if strcmp(nv.selection_ind, 'all') || isempty(nv.selection_ind)
    nv.selection_ind = 1:n_sessions; 
else
    assert(~isempty(nv.selection_ind))
    assert(isnumeric(nv.selection_ind))
    assert(all(floor(nv.selection_ind) == nv.selection_ind)) % Elements must be integer indices
end
n_selected_sessions = length(nv.selection_ind);

% Preallocate across selected sessions.
sessions_raw = cell(n_selected_sessions, 1);

% Load data from selected sessions. 
for i_session = 1:n_selected_sessions
    sessions_raw{i_session} ...
        = data_as_loaded.Raw.(all_session_ids{nv.selection_ind(i_session)});
end

% Get IDs of selected sessions.
selected_session_ids = all_session_ids(nv.selection_ind);


%% Check for NaNs in Lfold data

% Preallocate array whose i_th cell contains indices of regions that were
% removed for the i_th session.
nan_regions = cell(n_selected_sessions, 1);

% Preallocate array whose i_th cell contains the number of regions retained
% for the i_th session, after potentially dropping any regions.
n_regions_retained = NaN(n_selected_sessions, 1);

% Check all regions for each session, and remove regions with NaN values if
% desired. Optionally issue warning if regions removed.
for i_session = 1:n_selected_sessions
    nan_regions{i_session} = find(any(isnan(sessions_raw{i_session}.Lfold), 2));
    if nv.remove_nan && any(nan_regions{i_session})
        sessions_raw{i_session}.Lfold(nan_regions{i_session},:) = [];
        if nv.nan_warning
            warning('Region(s) %d from session %s have NaN values and will be removed.', ...
                    nan_regions{i_session}, selected_session_ids{i_session})
        end
    end
    
    % Note: in the future, as I become more familiar with these data, it is
    % likely that there are other arrays that will need to be similarly
    % modified.
    
    % Record number of regions retained for this session.
    n_regions_retained(i_session) = size(sessions_raw{i_session}.Lfold, 1);
end


%% Optionally plot Lfold from all selected sessions for current subject

if nv.plot      
    if nv.figure, figure; end

    % Plot all selected sessions for current subject.
    for i_session = 1:n_selected_sessions
    for i_region = 1:n_regions_retained(i_session)   
       subplot(n_selected_sessions, 1, i_session)
       hold on
       
       plot(i_region + sessions_raw{i_session}.Lfold(i_region,:)/nv.scale_factor, ...
            'linewidth', 2);   
       
       xlim(nv.window)
       xlabel('Time')
       ylabel('Region index')
       title(sprintf('Session %s', selected_session_ids{i_session}))     
    end
    end; clear i_session i_region

    % Resize figure window.
    set(gcf, 'Position', [100, 500, 1250, n_selected_sessions*400]);   
end

end

function activity = isolate_behavior(timeseries, session_raw, behavior)
% Pull out sections of Lfold, for each region, that correspond to the
% specified behavior of interest. 
% 
% PARAMETERS
% ----------
% timeseries  : m x n timeseries array, where m = number of variables, and
%               n = number of timepoints/bins.
% session_raw : Scalar struct output (possibly within cell array) from
%               extract_single_subj. Needed here are the fields .behavior
%               (list of all behaviors annotated within session) and
%               .Fstart and .Fstop (containing respectively, the start and
%               stop timepoints of the annotated behavior).
% behavior    : Behavior of interest which we would like to pull out.
%               String value, watch capitalization.
%
% RETURNS
% -------
% activity : m x p timeseries array, where the second array dim consists
%            of the concatenation of all m x p_i arrays, where i is in {1,
%            2, ... s}, with s = number of epochs with behavior of
%            interest, and m x p_i being the i_th timeseries array
%            annotated as behavior of interest.
%
% Author: Jonathan Chien. Refactored from local fx within xcorr_by_bhv on
% 5/31/22.


% Set offset to account for experimentalist reaction time. Note that this
% is in units of 40 ms. 
REACTION_TIME = 0;
session_raw.Fstart = session_raw.Fstart + REACTION_TIME;
session_raw.Fstop = session_raw.Fstop + REACTION_TIME;

% If any of the shifted timepoint indices now exceed the number of
% timepoints in timeseries, set to end of timeseries.
n_timepoints = size(session_raw.Lfold, 2);
session_raw.Fstart(session_raw.Fstart > n_timepoints) = n_timepoints;
session_raw.Fstop(session_raw.Fstop > n_timepoints) = n_timepoints;

% Ensure current session has requested behavior. If not throw exception
% with identifier in invoking function; it seems throwAsCaller will
% continue to pass the MException object up through calling functions until
% it reaches a try/catch block (or perhaps the highest level invoking
% function?).
if ~ismember(behavior, session_raw.behaviors)
    exception = MException("isolate_behavior:behavior_not_found", ...
                           "The specified behavior did not occur " + ...
                           "during the passed-in session.");
    throwAsCaller(exception);
end

% Get m x n matrix woi (windows of interest) where m = the number of epochs
% corresponding to the behavior of interest (the total number of epochs =
% length(session_raw.behaviors) and n = 2, where the 1st column is the
% start timepoint, and 2nd column is the stop timepoint.
all_windows = [session_raw.Fstart session_raw.Fstop];
woi = all_windows(ismember(session_raw.behaviors, behavior),:);

% Assemble an n_timepoints x n_regions matrix of neural activity, where the
% j_th column contains the concatenation of all windows of interest for the
% j_th region. Skip preallocation for increased readability.
activity = [];
for i_epoch = 1:size(woi, 1)
    activity = [activity, timeseries(:, woi(i_epoch,1):woi(i_epoch,2))];
end

end

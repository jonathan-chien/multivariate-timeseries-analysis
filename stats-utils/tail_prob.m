function p = tail_prob(observed,null,dim,nvp)
% Calculates the tail probability for an observation or set of observations
% under a null distribution or set of null distributions. Note that for the
% two-sided case, at least one source (An Introduction to Probability and
% Statistics, Second Edition, Rohatgi & Saleh) takes p values to be poorly
% defined for asymmetric null distributions, though this source states that
% doubling of the smaller tail is usually recommended by many authors in
% that case (also seems to be slightly more conservative than the absolute
% value method). Doubling may also perhaps be interpreted as a correction
% for two one-tailed tests.
%
% EXAMPLE
% -------
% size(observed) --> [20 5]
% size(null) --> [20 5 1000] (1000 draws from null distribution)
% p = tail_prob(observed, null, 3);
% size(p) --> [20 5]
%
% PARAMETERS
% ----------
% observed -- Array of dimensionality n whose elements are the observed
%             test statistics. E.g., a single test statistic corresponds to
%             n = 1, a vector of test statistics corresonds to n = 2, etc.,
%             though n may be any integer.
% null     -- Array of dimensionality n + 1 whose elements are draws from
%             the null distribution(s) of (each of) the passed in test
%             statistic(s) in "observed". Note that the array sizes and
%             orientations must be comptabile with array expansion, and all
%             array dimension sizes must match those of "observed" except
%             for the dim_th dimension (see dim below).
% dim      -- Integer scalar specifying the dimension along which are the
%             draws for each null distribution. If passed in empty, the
%             function will attempt to infer this value by comparing the
%             array dims of "observed" and "null". Note that this inference
%             will fail if the number of draws from the null
%             distribution(s) matches the size of any other array
%             dimension.
% Name-Value Pairs (nvp)
%   'type'       -- ("two-tailed" (default) | "right-tailed" |
%                   "left-tailed"), specify direction of test. Note that
%                   some authors consider the p value to be poorly defined
%                   for asymmetric null distributions, though doubling of
%                   the smaller tail is usually recommended in that case
%                   (which may also perhaps be viewed as a correction for
%                   testing two tails); this is the method adopted here.
%   'exact'      -- (1|0 default = 0), specifies whether to use calculation
%                   for an exact test or to correct for Monte Carlo
%                   simulation.
%   'correction' -- (1 | 0 (default)), specify whether or not (if true) to
%                   set any p values larger than 1, which may arise in some
%                   degenerate casee, to 1.
%
% RETURNS
% -------
% p -- Tail probability/ies of observation(s) under their respective null
%      distribution(s). Array size of "p" matches that of "observed", and
%      elements correspond one to one.
%
% Author: Jonathan Chien 8/20/21. Last edit: 6/2/22.

arguments
    observed
    null
    dim
    nvp.type = 'two-tailed'
    nvp.exact = false
    nvp.correction = false
end

% If user passed "dim" as an empty array, attempt to infer the correct
% value.
if isempty(dim)
    dim = setdiff(size(null), size(observed));
    if ~isscalar(dim)
        error("Unable to infer argument 'dim'. Please specify it.")
    end
end

% Calculate both right and left tails first. 
if ~nvp.exact 
    % Calculate both right and left tails with correction for
    % random permutations.
    nullDistrSize = size(null, dim);
    rightTail = squeeze(sum(observed <= null, dim) + 1) / (nullDistrSize + 1);
    leftTail = squeeze(sum(observed >= null, dim) + 1) / (nullDistrSize + 1);
            
elseif nvp.exact
    % Calculate both right and left tails without correction, if
    % exact permutations used.
    rightTail = mean(observed <= null, dim);
    leftTail = mean(observed >= null, dim);
    
else
    error("Invalid value for 'exact'.")
end

% Apply desired sidedness.
switch nvp.type
    case 'two-tailed'    
        % Double tails for all distributions and preallocate p.
        rightTailDoubled = rightTail * 2;
        leftTailDoubled = leftTail * 2;
        p = NaN(size(rightTail));
        
        % For each distribution, take smaller tail and double it. If tails
        % have equal mass, assign 1, as doubling in Monte Carlo case will
        % yield p > 1 due to correction.
        p(rightTail > leftTail) = leftTailDoubled(rightTail > leftTail);
        p(rightTail < leftTail) = rightTailDoubled(rightTail < leftTail);
        p(rightTail == leftTail) = ones(sum(rightTail==leftTail, 'all'), 1); 
        
    case 'right-tailed'
        p = rightTail;
        
    case 'left-tailed'
        p = leftTail;

    otherwise
        error("Invalid value for 'type'.")
end

% Set any p values greater than 1 to 1.
if nvp.correction, p(p > 1) = 1; end

end

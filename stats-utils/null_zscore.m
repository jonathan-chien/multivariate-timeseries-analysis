function [Z,nullMean,nullStdev] = null_zscore(observed,null,dim,varargin)
% Calculates the number of standard deviations from the mean of a null
% distribution that an observation/set of observations lie/lies.
%
% PARAMETERS
% ----------
% observed -- Array of dimensionality n whose elements are the observed
%             test statistics. E.g., a single test statistic corresponds to
%             n = 1, a vector of test statistics corresonds to n = 2, etc.
% null     -- Array of dimensionality n + 1 whose elements are draws from
%             the null distribution(s) of (each of) the passed in test
%             statistic(s) in observed. Note that the array sizes and
%             orientations must be comptabile with array expansion.
% dim      -- Specify array dimension of null along which the null values
%             for each variable lie. For example, if null is an m x n array
%             corresonding to n draws from each of m variable's associated
%             null distribution, dim should be set to 2.
% weight   -- Optional argument specifying weighting scheme for standard
%             deviation calculation.
% 
% RETURNS
% -------
% Z         -- Array of z-scores (number of standard deviations away from
%              the mean of the null distribution(s) that lie(s) the
%              observed value(s)). Array size of "z" matches that of
%              "observed", and elements correspond one to one.
% nullMean  -- Array of means of each null distribution. Shape matches that
%              of Z and observed.
% nullStdev -- Array of standard deviations of each null distribution.
%              Shape matches that of Z and observed.
%
% Author: Jonathan Chien 8/20/21. Last edit 6/9/22.


% Check if user passed in weight for standard deviation, else set to 0.
if ~isempty(varargin), weight = varargin{1}; else, weight = 0; end

% Compute mean and std of null distribution.
nullMean = squeeze(mean(null, dim));
nullStdev = squeeze(std(null, weight, dim));

% Ensure array sizes match. If nullMean matches observed, nullStdev should
% too.
assert(all(size(observed) == size(nullMean)), ...
       "After computing the mean of each null distribution, the resulting " + ...
       "array must match the shape of 'observed'.")

% Number of std from mean that observed lies.
Z = (observed - nullMean) ./ nullStdev;

end

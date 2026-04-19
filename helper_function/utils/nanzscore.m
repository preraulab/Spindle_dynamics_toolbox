%NANZSCORE  Z-score a numeric array while ignoring NaN entries
%
%   Usage:
%       [zscored, mu, sigma] = nanzscore(data, varargin)
%
%   Inputs:
%       data     : numeric - input array, may contain NaN -- required
%       varargin : additional arguments forwarded to zscore()
%
%   Outputs:
%       zscored : numeric - z-scored values of the non-NaN entries of data
%       mu      : double - mean of the non-NaN values
%       sigma   : double - standard deviation of the non-NaN values
%
%   Notes:
%       The returned zscored array has the same length as the non-NaN
%       subset of data, not the original input.
%
%   Example:
%       data = [1, 2, NaN, 4, 5];
%       [zscored, mu, sigma] = nanzscore(data);
%
%   See also: zscore, isnan
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

function [zscored, mu, sigma] = nanzscore(data, varargin)
    % Find non-NaN indices in the data
    inds = ~isnan(data);
    
    % Compute z-scores for non-NaN values
    [zscored, mu, sigma] = zscore(data(inds), varargin{:});
end



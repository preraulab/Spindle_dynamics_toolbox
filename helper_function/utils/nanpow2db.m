function ydB = nanpow2db(y)
%NANPOW2DB  Power-to-dB conversion with non-positive inputs mapped to NaN
%
%   Usage:
%       ydB = nanpow2db(y)
%
%   Inputs:
%       y : double - power value(s), any shape -- required
%
%   Outputs:
%       ydB : double - 10*log10(y), with NaN where y <= 0
%
%   Notes:
%       Adds and subtracts 300 to force integer results when y is an exact
%       power of 10. Adapted from The MathWorks pow2db with NaN-safe
%       handling of non-positive samples (MJP, 2020-02-07).
%
%   Example:
%       y1 = nanpow2db(2000/2);    % 30
%
%   See also: pow2db, db2pow, log10
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%ydB = 10*log10(y);
%ydB = db(y,'power');
% We want to guarantee that the result is an integer
% if y is a negative power of 10.  To do so, we force
% some rounding of precision by adding 300-300.

ydB = (10.*log10(y)+300)-300;
ydB(y(:)<=0) = nan;

% [EOF]
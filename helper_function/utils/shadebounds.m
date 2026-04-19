function [line_handle, bounds_handle] = shadebounds(x, mid, hi, lo, mid_color, shade_color, bounds_color, shade_alpha)
%SHADEBOUNDS  Plot a curve with shaded upper/lower bounds on a 2-D axis
%
%   Usage:
%       [line_handle, bounds_handle] = shadebounds(x, mid, hi, lo)
%       [line_handle, bounds_handle] = shadebounds(x, mid, hi, lo, mid_color, shade_color, bounds_color, shade_alpha)
%
%   Inputs:
%       x            : 1xN double - independent variable -- required
%       mid          : 1xN double - central curve -- required
%       hi           : 1xN double - upper bound of the shaded region -- required
%       lo           : 1xN double - lower bound of the shaded region -- required
%       mid_color    : color spec - color of the central curve (default: [0 0 0])
%       shade_color  : color spec - color of the shaded region (default: [.2 .2 .2])
%       bounds_color : color spec - color of the shaded region edges (default: [.9 .9 .9])
%       shade_alpha  : double - fill transparency (default: 0.1)
%
%   Outputs:
%       line_handle   : handle - handle to the central curve
%       bounds_handle : handle - handle to the fill object
%
%   Example:
%       x = 1:10; mid = x.^2; hi = mid + 6; lo = mid - 4;
%       shadebounds(x, mid, hi, lo);
%       shadebounds(x, mid, hi, lo, 'red', [1 0.8 0.8], 'black', 0.3);
%
%   See also: fill, plot
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

% Set defaults for optional input arguments
if nargin < 5
    mid_color = [0 0 0];              % Default curve color: black
end

if nargin < 6
    shade_color = [.2 .2 .2];      % Default shade color: dark gray
end

if nargin < 7
    bounds_color = [.9 .9 .9];        % Default bounds color: light gray
end

if nargin < 8
    shade_alpha = 0.1;                % Default shade transparency
end

% Make sure all vectors are pointing in the same direction
x = x(:)';
mid = mid(:)';
lo = lo(:)';
hi = hi(:)';

% Plot the curve and shaded bounds
hold on
bounds_handle = fill([x fliplr(x)], [lo fliplr(hi)], shade_color, 'EdgeColor', bounds_color, 'FaceAlpha', shade_alpha);
line_handle = plot(x, mid, 'LineWidth', 2, 'Color', mid_color);
end

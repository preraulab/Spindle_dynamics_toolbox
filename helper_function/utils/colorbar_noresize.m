function c = colorbar_noresize(varargin)
%COLORBAR_NORESIZE  Create a colorbar without resizing the parent axis
%
%   Usage:
%       c = colorbar_noresize
%       c = colorbar_noresize(ax, <colorbar arguments>)
%
%   Inputs:
%       ax : axes handle - parent axis for the colorbar (default: gca)
%
%   Outputs:
%       c : colorbar handle - handle to the created colorbar
%
%   Notes:
%       Preserves ax.Position across the call so the colorbar does not push
%       the axis. Any additional arguments are forwarded to colorbar().
%
%   See also: colorbar, axes
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

if nargin == 0 || ~isa(varargin{1},'matlab.graphics.axis.Axes')
    ax = gca;   % Use the current axis if no input argument is provided
else
    ax = varargin{1};
end

pos = ax.Position;  % Save the current position of the axis
c = colorbar(ax, varargin{2:end});  % Create the colorbar
ax.Position = pos;  % Restore the original position of the axis

end

function [ModelSpec] = specify_mdl(BinarySelect,InteractSelect,varargin)
%SPECIFY_MDL  Wrap user-provided model specifications into a ModelSpec struct
%
%   Usage:
%       ModelSpec = specify_mdl(BinarySelect, InteractSelect, 'Name', Value, ...)
%
%   Inputs:
%       BinarySelect   : 1x4 double - factor selector in fixed order [SOphase, stage, SOpower, history] -- required
%       InteractSelect : 1xn cell - interaction terms 'A:B' (case- and order-insensitive; : & - accepted) -- required
%
%   Name-Value Pairs:
%       'hist_choice'  : char - 'short' (15 s history, fast) or 'long' (90 s history, infraslow structure) (default: 'short')
%       'control_pt'   : 1xk double - spline control point locations (default: [0:15:90 120 150])
%       'binsize'      : double - point-process bin size in seconds (default: 0.1)
%       'hard_cutoffs' : 1x2 double - spindle frequency cutoffs in Hz (default: [12 16])
%
%   Outputs:
%       ModelSpec : struct - full model specification used throughout the toolbox
%
%   Example:
%       BinarySelect = [1, 1, 0, 1];
%       InteractSelect = {'stage:SOphase'};
%       ModelSpec = specify_mdl(BinarySelect, InteractSelect);
%
%   Notes:
%       Accompanies Chen et al., PNAS 2025. Valid interactions are
%       stage:SOphase, SOphase:SOpower, and stage:history.
%
%   See also: preprocessToDesignMatrix, build_design_mt, ismember_interaction
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;

% Required inputs
addRequired(p, 'BinarySelect', @(x) validateattributes(x, {'numeric', 'vector'}, {'real','nonempty','numel', 4}));
addRequired(p, 'InteractSelect', @(x) isempty(x) || iscell(x));

% Optional inputs
addOptional(p, 'hist_choice','short', @(x) isempty(x) || any(validatestring(x, {'short', 'long'})));
addOptional(p, 'control_pt', [0:15:90 120 150],@(x) isempty(x) || (isnumeric(x) && isvector(x) && all(isfinite(x)) && isreal(x)));
addOptional(p, 'binsize',0.1,@(x) validateattributes(x, {'numeric', 'vector'}, {'real','nonempty','nonnan','positive'}));
addOptional(p, 'hard_cutoffs',[12 16], @(x) validateattributes(x, {'numeric', 'vector'}, {'real','nonempty','nonnan','nonnegative'}));

parse(p, BinarySelect, InteractSelect,varargin{:});
parser_results = p.Results;

%% Save Model Spec

% Add factors selected into the model spec
AllFactors = {'SOphase','stage','SOpower','history'};
FactorSelect = AllFactors(parser_results.BinarySelect==1);      % Selected model components
ModelSpec = parser_results; 
ModelSpec.FactorSelect = strjoin(FactorSelect, ', ');

% Add hist_ord to the model spec
if ~isempty(parser_results.control_pt)
    ModelSpec.hist_ord = parser_results.control_pt(end);        % history order in bin
    ModelSpec.hist_ord_s = parser_results.control_pt(end)*parser_results.binsize; % hist order in sec
end

% Report selected components and interactions
disp(['Selected modeling components: ', strjoin(FactorSelect, ', ')]);
if isempty(parser_results.InteractSelect)
    disp('Selected interactions -> None')
elseif ~isempty(parser_results.InteractSelect)&&any(ismember_interaction(parser_results.InteractSelect,{'stage:SOphase','SOpower:SOphase','stage:history'}))
    disp(['Selected interactions -> ', strjoin(parser_results.InteractSelect, ', ')]);
else
    error('Interaction term is not properly specified. It must be specified as a cell array containing a subset of the valid options: {''stage:SOphase'', ''SOphase:SOpower'', ''stage:history''}');
end

end

  
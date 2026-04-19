function tf = ismember_interaction(targetInteraction, interactionsList)
%ISMEMBER_INTERACTION  Check whether interaction strings are members of a list, case- and order-insensitive
%
%   Usage:
%       tf = ismember_interaction(targetInteraction, interactionsList)
%
%   Inputs:
%       targetInteraction : char or cell of char - target interaction(s) of the form 'A:B' -- required
%       interactionsList  : cell of char - list of interactions to search within -- required
%
%   Outputs:
%       tf : 1xM logical - true where targetInteraction{m} matches any entry in interactionsList
%
%   Notes:
%       Comparisons are case-insensitive, order-insensitive (so 'A:B' and
%       'B:A' are equivalent), and treat ':', '&', and '-' as equivalent
%       separators.
%
%   Example:
%       interactions = {'stage:SOphase', 'stage:history'};
%       target = {'soPHASE&stage','stage:SOpower'};
%       tf = ismember_interaction(target, interactions);   % [1 0]
%
%   See also: ismember, specify_mdl, build_design_mt
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%%
    % Function to normalize the string
    normalize = @(s) strjoin(sort(split(regexprep(lower(s), '[:&-]', ':'), ':')), ':');
    
    % Normalize the list of interactions
    normalizedInteractionsList = cellfun(normalize, interactionsList, 'UniformOutput', false);
    
    % Normalize the target interaction(s)
    if iscell(targetInteraction)
        % If targetInteraction is a cell array, normalize each element
        normalizedTarget = cellfun(normalize, targetInteraction, 'UniformOutput', false);
    else
        % If targetInteraction is a single string, normalize directly
        normalizedTarget = {normalize(targetInteraction)};
    end
    
    % Check for matches for all normalized target interactions
    tf = ismember(normalizedTarget, normalizedInteractionsList);
end

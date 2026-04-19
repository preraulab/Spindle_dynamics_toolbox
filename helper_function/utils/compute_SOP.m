function [SOpower_norm, SOpower_times, SOpower_stages, norm_method, ptile] = compute_SOP(EEG, Fs, varargin)
%COMPUTE_SOP  Compute normalized slow-oscillation power from EEG
%
%   Usage:
%       [SOpower_norm, SOpower_times, SOpower_stages, norm_method, ptile] = ...
%           compute_SOP(EEG, Fs, 'Name', Value, ...)
%
%   Inputs:
%       EEG : 1xN double - time series EEG data -- required
%       Fs  : double - sampling frequency in Hz -- required
%
%   Name-Value Pairs:
%       'stage_vals'                : 1xS double - stage values (5=W,4=R,3=N1,2=N2,1=N3) (default: [])
%       'stage_times'               : 1xS double - stage onset times (default: [])
%       'SO_freqrange'              : 1x2 double - SO band in Hz (default: [0.3 1.5])
%       'SOpower_outlier_threshold' : double - z-score cutoff for SOpower outliers (default: 3)
%       'norm_method'               : char - 'pNshiftS', 'percent', 'proportion', or 'none' (default: 'p2shift1234')
%       'retain_Fs'                 : logical - upsample SOpower to EEG Fs (default: true)
%       'tapers'                    : 1x2 double - [TW, num_tapers] (default: [15 29])
%       'window_params'             : 1x2 double - [win_size, step_size] in seconds (default: [30 15])
%       'EEG_times'                 : 1xN double - EEG sample times (default: (0:N-1)/Fs)
%       'time_range'                : 1x2 double - restrict normalization to this window (default: [EEG_times(1), end])
%       'isexcluded'                : 1xN logical - artifact mask (default: all false)
%
%   Outputs:
%       SOpower_norm   : 1xM double - normalized SOpower time series (row vector, possibly upsampled)
%       SOpower_times  : 1xM double - times of the SOpower samples (s)
%       SOpower_stages : 1xM double - stage value at each SOpower sample (or true if staging not provided)
%       norm_method    : char - resolved normalization method name
%       ptile          : double or 1x2 double - percentile value(s) used for normalization (empty for 'none'/'proportion')
%
%   Notes:
%       pNshiftS parses to percentile N computed across the union of stages
%       in S (digits). Artifact samples are replaced with NaN before the
%       multitaper estimate via compute_mtspect_power.
%
%   See also: compute_mtspect_power, nanzscore, multitaper_spectrogram_mex
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Parse input
%Input Error handling
p = inputParser;

%Stage info
addOptional(p, 'stage_vals', [], @(x) validateattributes(x, {'double', 'single'}, {'real'}));
addOptional(p, 'stage_times', [], @(x) validateattributes(x, {'numeric', 'vector'}, {'real'}));

%SOpower settings
addOptional(p, 'SO_freqrange', [0.3, 1.5], @(x) validateattributes(x, {'numeric', 'vector'}, {'real', 'finite', 'nonnan'}));
addOptional(p, 'SOpower_outlier_threshold', 3, @(x) validateattributes(x,{'numeric'}, {'scalar'}));
addOptional(p, 'norm_method', 'p2shift1234', @(x) validateattributes(x, {'char', 'numeric'},{}));
addOptional(p, 'retain_Fs', true, @(x) validateattributes(x,{'logical'},{}));
addOptional(p, 'tapers', [15 29], @(x) validateattributes(x,{'numeric', 'vector'}, {'numel',2}));
addOptional(p, 'window_params', [30 15], @(x) validateattributes(x,{'numeric', 'vector'}, {'numel',2}));

%EEG time settings
addOptional(p, 'EEG_times', [], @(x) validateattributes(x, {'numeric', 'vector'},{'real','finite','nonnan'}));
addOptional(p, 'time_range', [], @(x) validateattributes(x, {'numeric', 'vector'},{'real','finite','nonnan'}));
addOptional(p, 'isexcluded', [], @(x) validateattributes(x, {'logical', 'vector'},{}));

parse(p,varargin{:});
parser_results = struct2cell(p.Results); %#ok<NASGU>
field_names = fieldnames(p.Results);

eval(['[', sprintf('%s ', field_names{:}), '] = deal(parser_results{:});']);

if isempty(EEG_times) %#ok<*NODEF>
    EEG_times = (0:length(EEG)-1)/Fs;
else
    assert(length(EEG_times) == size(EEG,2), 'EEG_times must be the same length as EEG');
end

if isempty(time_range)
    time_range = [min(EEG_times), max(EEG_times)];
else
    assert( (time_range(1) >= min(EEG_times)) & (time_range(2) <= max(EEG_times)), 'time_range cannot be outside of the time range described by "EEG_times"');
end

if isempty(isexcluded)
    isexcluded = false(size(EEG,2),1);
else
    assert(length(isexcluded) == size(EEG,2),'isexcluded must be the same length as EEG');
end

%% Compute SO power
% Replace artifact timepoints with NaNs
nanEEG = EEG;
nanEEG(isexcluded) = nan;

% Now compute SOpower
[SOpower, SOpower_times] = compute_mtspect_power(nanEEG, Fs, 'freq_range', SO_freqrange,'tapers', tapers, 'window_params', window_params);

SOpower_times = SOpower_times + EEG_times(1); % adjust the time axis to EEG_times

% Compute SOpower stage
if ~isempty(stage_vals) && ~isempty(stage_times)
    SOpower_stages = interp1(stage_times, stage_vals, SOpower_times, 'previous');
else
    SOpower_stages = true;
end

% Exclude outlier SOpower that usually reflect isexcluded
SOpower(abs(nanzscore(SOpower)) >= SOpower_outlier_threshold) = nan;


%% Normalize SO power

%Handle shift inputs
if strcmpi(norm_method,'shift')
    shift_ptile = 2;
    shift_stages = 1:4;
elseif regexp(norm_method,'p*+shift*')
    shift_ptile = str2double(norm_method(2:strfind(norm_method,'shift')-1));
    shift_stages = unique((norm_method(strfind(norm_method,'shift')+5:end)) - '0');

    if isempty(shift_stages)
        shift_stages = 1:4;
    end

    norm_method = 'shift';
end

%Check shift
assert(shift_ptile >= 0 && shift_ptile <= 100, 'Shift percentile must be between 0 and 100');

switch norm_method
    case {'proportion', 'normalized'}
        [proppower, ~] = computeMTSpectPower(nanEEG, Fs, 'freq_range', [freq_range(1), freq_range(2)]);
        SOpower_norm = db2pow(SOpower)./db2pow(proppower);
        ptile = [];

    case {'percentile', 'percent', '%', '%SOP'}
        low_val =  1;
        high_val =  99;
        ptile = prctile(SOpower(SOpower_times>=time_range(1) & SOpower_times<=time_range(2)), [low_val, high_val]);
        SOpower_norm = SOpower-ptile(1);
        SOpower_norm = SOpower_norm/(ptile(2) - ptile(1));  % Note: is this computation correct? 

    case {'shift'}
        if islogical(SOpower_stages) && SOpower_stages
            SOpower_stages_valid = true(size(SOpower_stages));
        else
            SOpower_stages_valid = ismember(SOpower_stages, shift_stages);
        end
        ptile = prctile(SOpower(SOpower_times>=time_range(1) & SOpower_times<=time_range(2) & SOpower_stages_valid), shift_ptile);
        SOpower_norm = SOpower-ptile(1);

    case {'absolute', 'none'}
        SOpower_norm = SOpower;
        ptile = [];

    otherwise
        error(['Normalization method "', norm_method, '" not recognized']);
end

% Make the output SOpower_norm a row vector 
SOpower_norm = SOpower_norm';

%% (Optional) Upsample to EEG sampling rate 
if retain_Fs
    SOpower_norm_notnan = SOpower_norm(~isnan(SOpower_norm));
    SOpower_norm = interp1([EEG_times(1), SOpower_times(~isnan(SOpower_norm)), EEG_times(end)],...
       [SOpower_norm_notnan(1), SOpower_norm_notnan, SOpower_norm_notnan(end)], EEG_times);
    SOpower_norm(isexcluded) = nan;
    SOpower_times = EEG_times;
    if ~isempty(stage_vals) && ~isempty(stage_times)
        SOpower_stages = interp1(stage_times, stage_vals, SOpower_times, 'previous');
    else
        SOpower_stages = true;
    end
end

end

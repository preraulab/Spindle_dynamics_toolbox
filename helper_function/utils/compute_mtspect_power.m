function [SO_power, stimes, sfreqs] = compute_mtspect_power(varargin)
%COMPUTE_MTSPECT_POWER  Compute band power from the multitaper spectrogram of a time series
%
%   Usage:
%       [SO_power, stimes, sfreqs] = compute_mtspect_power(data, Fs, 'Name', Value, ...)
%
%   Inputs:
%       data : 1xN double - time series -- required
%       Fs   : double - sampling frequency in Hz -- required
%
%   Name-Value Pairs:
%       'freq_range'       : 1x2 double - band to integrate in Hz (default: [0.3 1.5])
%       'tapers'           : 1x2 double - [time-halfbandwidth product, number of tapers] (default: [15 29])
%       'window_params'    : 1x2 double - [window size (s), step size (s)] (default: [30 15])
%       'smoothing_method' : char - smoothdata method or 'none' (default: 'none')
%       'smoothing_param'  : double - smoothing window in seconds (default: 300)
%       'interp_times'     : 1xM double - times at which to interpolate the output (default: [])
%       'verbose'          : logical - print diagnostics (default: false)
%
%   Outputs:
%       SO_power : 1xT double - band-integrated power in dB (or interp_times when specified)
%       stimes   : 1xT double - center times of the spectrogram bins (s)
%       sfreqs   : 1xF double - spectrogram frequency bins (Hz)
%
%   Notes:
%       Uses multitaper_spectrogram_mex internally; SO_power is converted to
%       dB via nanpow2db with NaNs preserved across smoothing gaps.
%
%   See also: multitaper_spectrogram_mex, nanpow2db, smoothdata
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Parse input
%Input Error handling
p = inputParser;

addRequired(p,'data',@(x) validateattributes(x,{'numeric', 'vector'},{'nonempty'}));
addRequired(p,'Fs',@(x) validateattributes(x,{'numeric'},{'nonempty','numel',1}));
addOptional(p,'freq_range',[.3 1.5], @(x) validateattributes(x,{'numeric', 'vector'},{'numel',2}));
addOptional(p,'tapers',[15 29],@(x) validateattributes(x,{'numeric', 'vector'},{'numel',2}));
addOptional(p,'window_params',[30 15], @(x) validateattributes(x,{'numeric', 'vector'},{'numel',2}));
addOptional(p,'smoothing_method','none', @(x) validateattributes(x,{'char'},{}));
addOptional(p,'smoothing_param', 60*5, @(x) validateattributes(x,{'numeric'},{'numel',1}));
addOptional(p,'interp_times',[], @(x) validateattributes(x,{'numeric', 'vector'},{'real'}));
addOptional(p,'verbose',false, @(x) validateattributes(x,{'logical'},{}));
parse(p,varargin{:});
parser_results = struct2cell(p.Results);
field_names = fieldnames(p.Results);

eval(['[', sprintf('%s ', field_names{:}), '] = deal(parser_results{:});']);

%% Compute SO-power
%Compute power using the MTS (data, Fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose)

[SO_spect, stimes, sfreqs] = multitaper_spectrogram_mex(data, Fs, freq_range, tapers, window_params, [], [], [], false, verbose);

%Compute dt
dt = stimes(2)-stimes(1);
df = sfreqs(2) - sfreqs(1);

%Takes the total power and converts to dB
SO_power = nanpow2db(sum(SO_spect,1)*df)';

%% Smooth data
if ~strcmpi(smoothing_method, 'none') && ~isempty(smoothing_param) && smoothing_param>0
    
    if verbose
        disp(['Smoothing using ' smoothing_method ' with parameter ' num2str(smoothing_param)]);
    end
    
    %Get bad indices
    bad_inds = ~isfinite(SO_power);
    
    %Interpolate big gaps in data
    t = 1:length(SO_power);
    data_fixed = interp1([0, t(~bad_inds), length(SO_power)+1], [0; SO_power(~bad_inds); 0], t)';
    
    smooth_samples   = smoothing_param/dt;        %Time in samples
    
    SO_power = smoothdata(data_fixed, smoothing_method, smooth_samples, 'omitnan' );
    
    %Return the bad values
    SO_power(bad_inds) = nan;
end

%% Interpolate data
if ~isempty(interp_times)
    SO_power = interp1(stimes, SO_power, interp_times);
end
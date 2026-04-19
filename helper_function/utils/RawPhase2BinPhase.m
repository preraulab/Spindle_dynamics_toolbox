function [SO_phase] = RawPhase2BinPhase(phase_raw,Fs,spindletime_raw,spindletrain,domaindivide)
%RAWPHASE2BINPHASE  Align raw SO phase to a binned spindle train (exact phase on events, bin mean elsewhere)
%
%   Usage:
%       SO_phase = RawPhase2BinPhase(phase_raw, Fs, spindletime_raw, spindletrain, domaindivide)
%
%   Inputs:
%       phase_raw       : 1xM double - raw SO phase time series (rad, will be unwrapped internally) -- required
%       Fs              : double - sampling frequency of phase_raw in Hz -- required
%       spindletime_raw : 1xK double - raw spindle event times (s) -- required
%       spindletrain    : 1xB double - binned spindle indicator train -- required
%       domaindivide    : 1x(B+1) double - bin edges in seconds -- required
%
%   Outputs:
%       SO_phase : Bx1 double - binned SO phase; exact raw phase when a spindle occurs in a bin,
%                  otherwise the circular mean of phase_raw over that bin
%
%   See also: rawToBinData, eegToEventSignal, unwrap
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% SO phase   

phase_raw = unwrap(phase_raw);                          % unwrap the phase for taking avg
phase_time_raw = 1/Fs:1/Fs:length(phase_raw)/Fs;        % phase time

% Construct SO phase train
% 1) For phase train with spindle events, find the exact phase for each spindle 
spindle_idx = find(histcounts(spindletime_raw,phase_time_raw));
spindle_phase = phase_raw(spindle_idx);

% 2) For phase train without spindle events, take the mean of those falling into each bin 
bincount_phase = histcounts(phase_time_raw,domaindivide)';  
bincount_phase2 = cumsum(bincount_phase);
bincount_phase3 = bincount_phase2 - (bincount_phase-1);

phase_train1 = zeros(length(spindletrain),1);
for i = 1:length(phase_train1)
    phase_train1(i) = mean(phase_raw(bincount_phase3(i):bincount_phase2(i)));
end

% 3) Replace the exact phase when spindle occurs
SO_phase = phase_train1;
total_spikes = length(find(spindletrain==1));
SO_phase(spindletrain==1) = spindle_phase(1:total_spikes);

end
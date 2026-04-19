function [BinData] = rawToBinData(res_table,Fs,binsize,hard_cutoffs,hist_choice)
%RAWTOBINDATA  Convert extracted event/signal data into a binned struct ready for GLM fitting
%
%   Usage:
%       BinData = rawToBinData(res_table, Fs, binsize, hard_cutoffs, hist_choice)
%
%   Inputs:
%       res_table    : table - extracted spindle event info and SO signals -- required
%       Fs           : double - sampling frequency in Hz -- required
%       binsize      : double - point-process bin size in seconds -- required
%       hard_cutoffs : 1x2 double - frequency range (Hz) for fast spindle selection -- required
%       hist_choice  : char - 'short' (15 s) or 'long' (90 s, infraslow) history window -- required
%
%   Outputs:
%       BinData : struct - binned model inputs (y, sta, sop, phase, sp_hist, sp, RW,
%                 real_time, spindletime_raw, isis, N2_rate, N3_rate, Fs)
%
%   Notes:
%       The recording is trimmed to a 5-minute (300 s) buffer around the
%       first and last NonWake epoch. History uses cardinal spline control
%       points [0:15:90 120 hist_ord] (short) or
%       [0:15:90 120 150:100:750 hist_ord] (long). Accompanies Chen et al.,
%       PNAS 2025.
%
%   See also: preprocessToDesignMatrix, ModifiedCardinalSpline, Hist, RawPhase2BinPhase, stageraw2bin
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Extract fast spindle info
fast_inds = (res_table.peak_freqs{:} >= hard_cutoffs(1)) & (res_table.peak_freqs{:}< hard_cutoffs(2));
eventtime_all = res_table.peak_ctimes{:};             % All peak times

spindletime_raw = eventtime_all(fast_inds);           % Fast spindle times
if iscell(res_table.t)
    stagetime_raw = cell2mat(res_table.t);              
    stages_raw = cell2mat(res_table.stages);  
else 
    stagetime_raw = res_table.t;                          % Stage time
    stages_raw = res_table.stages;                        % Stage values
end
phase_raw = cell2mat(res_table.SOphase);              % SO phase 
sop_raw = cell2mat(res_table.SOpow);                  % SO power  

%% Convert data to binned data 
epoch = stagetime_raw(2)-stagetime_raw(1);                    % length of each epoch
totaltime = stagetime_raw(end) - stagetime_raw(1) + epoch;  
timestamps = 0:binsize:totaltime;    

% Response: binary spindle train
spindletrain= histcounts(spindletime_raw,timestamps)';       % distribute event times to intervals

% Predictors
% 1) sleep stage: stagetrain, stage matrix 
[stage,stagetrain] = stageraw2bin(stages_raw,stagetime_raw,binsize,epoch,timestamps);

% 2) SO phase   
phase_time_raw = 1/Fs:1/Fs:length(phase_raw)/Fs;        % phase time
SO_phase = RawPhase2BinPhase(phase_raw,Fs,spindletime_raw,spindletrain,timestamps);% SO phase raw in bin

% 3) SO power
SO_power = interp1(phase_time_raw, sop_raw, binsize:binsize:totaltime)';  % SO power raw in bin    

%% Cut data by the 5-min (300 sec) lights on and off
buffer_length = 300;     % 300 sec, 5 min               

% lights off idx: max(5 min before NonWake stage,1)
light_off_idx = max(find(stagetrain~=5,1,'first')-buffer_length/binsize,1);
% lights on idx: min(5 min after NonWake stage,length(spiketrain))
light_on_idx = min(find(stagetrain~=5,1,'last') + buffer_length/binsize, length(stagetrain)); 

spindletrain_cut = spindletrain(light_off_idx:light_on_idx);
stage_cut = stage(light_off_idx:light_on_idx,:);
stagetrain_cut = stagetrain(light_off_idx:light_on_idx);
sop_cut = SO_power(light_off_idx:light_on_idx);
phase_cut = SO_phase(light_off_idx:light_on_idx);
timestamps_cut = timestamps(light_off_idx:light_on_idx);

%% Set history order
switch lower(hist_choice)
    case {'short'}
          hist_ord = 15/binsize;   % 15 sec history
          control_pt = [0:15:90 120 hist_ord]; 
    case {'long'}
          hist_ord = 90/binsize;   % 90 sec history
          control_pt = [0:15:90 120 150:100:750 hist_ord]; 
end

%% Build up history component
indi_hist = Hist(hist_ord,spindletrain_cut);  % indicator-based history 

% Cardinal spline based history
s = .5;                                               % spline tension parameter          
[sp] = ModifiedCardinalSpline(hist_ord,control_pt,s);  % spline coefficients
sp_hist = indi_hist*sp;                               % spline-based history

%% Construct design matrix and response aligned with hist component 

y = spindletrain_cut(hist_ord+1:end);
sta = stage_cut(hist_ord+1:end,:);
sop = sop_cut(hist_ord+1:end,:);
phase = phase_cut(hist_ord+1:end,:);
real_time = timestamps_cut(hist_ord+1:end)';
RW = sta(:,4)+sta(:,5); % REM/Wake indicator

% spindle density in mins
N2_rate = (length(find(spindletrain_cut>=1&stagetrain_cut==2))/length(spindletrain_cut(stagetrain_cut==2)))*(60/binsize);
N3_rate = (length(find(spindletrain_cut>=1&stagetrain_cut==1))/length(spindletrain_cut(stagetrain_cut==1)))*(60/binsize);


%% Save info and output a struct

BinData = struct();
BinData.y = y;
BinData.sta = sta;
BinData.sop = sop;
BinData.phase = phase;
BinData.sp_hist = sp_hist;
BinData.sp = sp;
BinData.RW = RW;
BinData.real_time = real_time;
BinData.spindletime_raw = spindletime_raw;
BinData.isis = diff(spindletime_raw);
BinData.N2_rate = N2_rate;
BinData.N3_rate = N3_rate;
BinData.Fs = Fs;
end

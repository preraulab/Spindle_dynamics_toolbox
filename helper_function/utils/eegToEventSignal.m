function [res_table]  = eegToEventSignal(EEG,Fs,stage_val,stage_time)
%EEGTOEVENTSIGNAL  Preprocess EEG to extract spindle events and slow-oscillation power/phase
%
%   Usage:
%       res_table = eegToEventSignal(EEG, Fs, stage_val, stage_time)
%
%   Inputs:
%       EEG        : 1xN double - raw EEG data -- required
%       Fs         : double - sampling frequency in Hz -- required
%       stage_val  : 1xS double - stage values (1=N3, 2=N2, 3=N1, 4=REM, 5=Wake) -- required
%       stage_time : 1xS double - stage onset times (s) -- required
%
%   Outputs:
%       res_table : table - event info and signals with key columns:
%           peak_ctimes : detected event central times (s)
%           peak_freqs  : detected event frequency (Hz)
%           SOpow       : slow oscillation power (normalized)
%           SOphase     : slow oscillation phase (rad, unwrapped on input)
%
%   Notes:
%       Pipeline: detect_artifacts -> TF_peak_detection -> compute_SOP
%       (p2shift1234 norm) -> SO band-pass filter + Hilbert for phase.
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: detect_artifacts, TF_peak_detection, compute_SOP, rawToBinData
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Presettings

N = length(EEG);
t = (0:N-1)/Fs;  % time stamp
SOpower_outlier_threshold = 3;
norm_method = 'p2shift1234';
stage_val_t = interp1(stage_time,stage_val, t, 'previous');
chan_names = "central"; % if chan_name is not provided
night = 1; % if night is not provided
ID = 1; % if subj ID is not provided

% check EEG polarity (e.g., MESA has polarity issue, needs to be flipped) 
% EEG = -EEG;

%% Precompute SO bandpass filter 
SO_freqrange = [0.4, 1.5];
d = designfilt('bandpassiir', ...       % Response type
       'StopbandFrequency1',SO_freqrange(1)-0.1, ...    % Frequency constraints
       'PassbandFrequency1',SO_freqrange(1), ...
       'PassbandFrequency2',SO_freqrange(2), ...
       'StopbandFrequency2',SO_freqrange(2)+0.1, ...
       'StopbandAttenuation1',60, ...   % Magnitude constraints
       'PassbandRipple',1, ...
       'StopbandAttenuation2',60, ...
       'DesignMethod','ellip', ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',Fs);

%% Run spindle detection and compute SO
% Artifact detection
artifacts = detect_artifacts(EEG, Fs);

% % TF peak detection
[spindle_table,~] = TF_peak_detection(EEG,Fs,[t; stage_val_t]','artifact_detect',artifacts);

% Normalized SO Power 
[SOP_norm0, SOP_time0,~,~] = compute_SOP(EEG, Fs, 'stage_vals', stage_val_t, 'stage_times', t,...
            'SO_freqrange', SO_freqrange, 'SOpower_outlier_threshold', SOpower_outlier_threshold, 'norm_method', norm_method,...
            'EEG_times', t, 'time_range', [t(1), t(end)], 'isexcluded', artifacts,'window_params', [4 1],'tapers', [2 3]);

[SOP_norm,~] = fillmissing(SOP_norm0,'linear','SamplePoints',SOP_time0);

% SO filtered signal, SO phase and amplitude
 SO_filt = filtfilt(d,EEG);         % filtered SO
 hilbert_sig = hilbert(SO_filt);    % hilbert signal
 phase_raw = angle(hilbert_sig);    % SO phase   
% amp_raw = abs(hilbert_sig);       % SO amplitude

%% Append data
% Intialize res_table to store data
varTypes = ["double","double","string","cell","cell",...
            "cell","cell","cell","cell","cell",...
            "cell","cell","cell","cell","cell","cell",...
            "cell","cell"];
varNames = ["ID","night","channel","stages","t", ...
           "peak_stimes","peak_etimes","peak_ctimes","peak_maxprom_times","peak_stages",...
           "peak_freqs","peak_proms","peak_durs","peak_freqHi","peak_freqLo","peak_vol",...
           "SOpow","SOphase"];
sizetable = [1,length(varNames)];
res_table = table('Size',sizetable,'VariableNames',varNames,'VariableTypes',varTypes);

% Store data
res_table.ID = ID;
res_table.night = night;
res_table.channel = chan_names;
res_table.stages{:} = stage_val_t;
res_table.t{:} = t;

res_table.peak_ctimes{:} = (spindle_table.Start_Time + spindle_table.End_Time)/2 ; 
res_table.peak_stimes{:} = spindle_table.Start_Time ;
res_table.peak_etimes{:} = spindle_table.End_Time ;
res_table.peak_maxprom_times{:} = spindle_table.Max_Prominence_Time;
res_table.peak_freqs{:} = spindle_table.Freq_Central;
res_table.peak_stages{:} = spindle_table.Stage;
res_table.peak_proms{:} = spindle_table.Prominence_Value;
res_table.peak_durs{:} = spindle_table.Duration;
res_table.peak_freqHi{:} = spindle_table.Freq_High;
res_table.peak_freqLo{:} = spindle_table.Freq_Low;
res_table.peak_vol{:} = spindle_table.Peak_Volume;


res_table.SOpow{:} = SOP_norm;
res_table.SOphase{:} = phase_raw;

% res_table.SOpow0{:} = SOP_norm0;
% res_table.SOamp{:} = amp_raw;
% res_table.SOfilt{:} = SO_filt;

%save(fullfile(savefolder_name, ['Example',num2str(ID),'_night',num2str(night),'_EventSignal.mat']), 'res_table', '-v7.3');
%save(['Example',num2str(ID),'_nt',num2str(night),'_EventSignal.mat'], 'res_table', '-v7.3');

end
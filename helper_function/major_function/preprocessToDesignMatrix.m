function [X, BinData, res_table] = preprocessToDesignMatrix(EEG, Fs, stage_val, stage_time,ModelSpec)
%PREPROCESSTODESIGNMATRIX  Preprocess raw EEG data into a design matrix and binned data
%
%   Usage:
%       [X, BinData, res_table] = preprocessToDesignMatrix(EEG, Fs, stage_val, stage_time, ModelSpec)
%
%   Inputs:
%       EEG        : 1xN double - raw EEG time series -- required
%       Fs         : double - sampling frequency in Hz -- required
%       stage_val  : 1xS double - stage values (1,2,3,4,5 = N3,N2,N1,REM,Wake) -- required
%       stage_time : 1xS double - stage onset times (s) -- required
%       ModelSpec  : struct - model specifications from specify_mdl -- required
%
%   Outputs:
%       X         : TxP double - design matrix (size depends on data length and ModelSpec)
%       BinData   : struct - binned data aligned to binsize
%       res_table : table - extracted event info and signals (peak_ctimes, peak_freqs, SOpow, SOphase, ...)
%
%   Notes:
%       Pipeline: eegToEventSignal -> rawToBinData -> build_design_mt.
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: eegToEventSignal, rawToBinData, build_design_mt, specify_mdl
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Extract spindle events and their properties

binsize = ModelSpec.binsize;
hard_cutoffs = ModelSpec.hard_cutoffs; 
hist_choice = ModelSpec.hist_choice;

disp('May take minutes ... Extracting spindles ...');
[res_table] = eegToEventSignal(EEG, Fs, stage_val, stage_time);

%% Convert event data to binned data
[BinData] = rawToBinData(res_table, Fs, binsize, hard_cutoffs, hist_choice);

%% Build design matrix for the specified model
X = build_design_mt(ModelSpec, BinData);
    
end

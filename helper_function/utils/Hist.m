function [SpikeHist] = Hist(lag,spiketrain)
%HIST  Build an indicator-basis history matrix from a binary spike train
%
%   Usage:
%       SpikeHist = Hist(lag, spiketrain)
%
%   Inputs:
%       lag        : integer - history lag in bins -- required
%       spiketrain : 1xN double/logical - binary event train -- required
%
%   Outputs:
%       SpikeHist : (N-lag)xlag double - indicator-basis history matrix;
%                   column i holds spiketrain shifted by i bins
%
%   See also: ModifiedCardinalSpline, build_design_mt
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

SpikeHist = zeros(length(spiketrain)-lag,lag);

for i = 1:lag
    SpikeHist(:,i) = spiketrain(lag-i+1:end-i);     % Shift the train
end

end
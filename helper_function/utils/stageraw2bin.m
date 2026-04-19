function [stage,stagetrain] = stageraw2bin(stages_raw,stagetime_raw,bin,epoch,domaindivide)
%STAGERAW2BIN  Convert (stage_time, stage_val) hypnogram into a binned stage train and one-hot stage matrix
%
%   Usage:
%       [stage, stagetrain] = stageraw2bin(stages_raw, stagetime_raw, bin, epoch, domaindivide)
%
%   Inputs:
%       stages_raw    : 1xS double - raw stage values (1=N3, 2=N2, 3=N1, 4=REM, 5=Wake) -- required
%       stagetime_raw : 1xS double - raw stage onset times (s) -- required
%       bin           : double - target bin size (s) -- required
%       epoch         : double - stage epoch length (s) -- required
%       domaindivide  : 1x(B+1) double - bin edges -- required
%
%   Outputs:
%       stage      : Bx5 double - one-hot stage matrix in [N1 N2 N3 REM Wake] column order
%       stagetrain : Bx1 double - per-bin stage code
%
%   Notes:
%       When bin <= epoch, samples are expanded by the epoch/bin ratio;
%       otherwise stages are aggregated per-bin via histcounts.
%
%   See also: rawToBinData, histcounts
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

% sleep stage train
if bin<= epoch
    stagetrain = zeros(length(domaindivide)-1,1);
   for i = 1:length(stagetime_raw)
      stagetrain(epoch/bin*(i-1)+1:(epoch/bin)*i) = stages_raw(i); 
   end
else 
    stagetrain = stages_raw(histcounts(domaindivide,stagetime_raw)==1)';   % sleep stage train complete
end

% Transfer stage train to stage matrix
stage = zeros(length(domaindivide)-1,5);                       
stage(stagetrain==3,1)=1;                                   % N1
stage(stagetrain==2,2)=1;                                   % N2
stage(stagetrain==1,3)=1;                                   % N3
stage(stagetrain==4,4)=1;                                   % REM
stage(stagetrain==5|stagetrain==0|isnan(stagetrain),5)=1;   % WAKE

end
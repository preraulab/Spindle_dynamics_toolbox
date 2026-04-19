function [stage_prop,timestamps] = compute_stage_prob(stagetrain,sop,sop_bin_divide)
%COMPUTE_STAGE_PROB  Compute per-bin stage proportions for a feature (e.g., SOpower)
%
%   Usage:
%       [stage_prop, timestamps] = compute_stage_prob(stagetrain, sop, sop_bin_divide)
%
%   Inputs:
%       stagetrain     : 1xN double - interpolated hypnogram at the resolution of sop -- required
%       sop            : 1xN double - feature value per sample (e.g., SOpower) -- required
%       sop_bin_divide : 1xB double - bin edges along the feature axis -- required
%
%   Outputs:
%       stage_prop : Bx5 double - per-bin stage proportions [N3, N2, N1, REM, Wake]
%       timestamps : 1xB double - bin centers (useful for overlaying with the feature axis)
%
%   Notes:
%       Each bin collects samples whose sop falls strictly between
%       consecutive edges; the histogram over stage values 1..5 is
%       normalized to 1 per bin.
%
%   See also: compute_SOP, histcounts
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%Compute the stage probability and components
bin_size = sop_bin_divide(2) - sop_bin_divide(1);
num_bins = length(sop_bin_divide);
stage_prop = zeros(num_bins, 5);

%Loop through each bin and get the associated stages
for ii = 1:num_bins-1
    inds = sop>sop_bin_divide(ii)&sop<sop_bin_divide(ii+1);
    hist_vals = histcounts(stagetrain(inds),.5:1:5.5)'; %Create histogram bin edges that include 1:5
    stage_prop(ii,:) = hist_vals/sum(hist_vals); % N3,N2,N1,REM,Wake
end

% The center of bin can be the timestamps to plot with SOP
timestamps = sop_bin_divide + 0.5*bin_size;


end
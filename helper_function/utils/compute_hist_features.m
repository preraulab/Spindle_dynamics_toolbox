function [hist_features] = compute_hist_features(xlag,yhat,yl,yu)
%COMPUTE_HIST_FEATURES  Compute refractory, excited, and peak features from a history-modulation curve
%
%   Usage:
%       hist_features = compute_hist_features(xlag, yhat, yl, yu)
%
%   Inputs:
%       xlag : nx1 double - history time lag (s) -- required
%       yhat : nx1 double - history modulation curve (rate multiplier) -- required
%       yl   : nx1 double - 95 percent CI lower bound -- required
%       yu   : nx1 double - 95 percent CI upper bound -- required
%
%   Outputs:
%       hist_features : struct with fields
%           ref_period : refractory period (s) - first time the upper bound exceeds 1 with yhat still < 1
%           exc_period : excited period (s) - length of the run of yl >= 1 that contains the first valid peak
%           p_time     : peak time (s)
%           p_height   : peak height (clamped to <= 1 if no significant peak)
%           AUC_is     : rate multiplier averaged over the 40-70 s infraslow window (only when xlag spans > 1500 bins)
%
%   See also: plot_hist_curve, consecutive_runs, findpeaks
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Compute hist features 
xbin = xlag(2)-xlag(1);

% compute upper bound gradient
yu_g = gradient(yu);

% sig ref period
idx_over1= find(yu>1&yu_g>=0&yhat<1,1);
ref_period = xlag(idx_over1);

% find all peaks
[p_h,p_locs,~,~] = findpeaks(yhat,'MinPeakHeight',1,'WidthReference','halfheight');

%---- choose the first sig peak----
% 1) compute periods with LB >= 1
[run_len, run_inds, yes_vector] = consecutive_runs(yl>=1, 1, inf, 1);

% 2) valid peaks
valid_peak = yes_vector(p_locs)== 1;

% 3) Find first valid peak idx
valid_peak_idx = find(valid_peak==1,1);

% 4) if LB is always< 1, choose first peak as the peak, set pw = 0, compute pt and ph
if isempty(valid_peak_idx)||xlag(p_locs(valid_peak_idx))>=8  
    [m_h,m_idx] = max(yhat(1:800));
    p_height = min(1,m_h);
    p_time = xlag(m_idx);
    exc_period = 0;

else % Find corresponding runs contains valid peak, use run length as peak width
    p_height = p_h(valid_peak_idx);
    p_time = xlag(p_locs(valid_peak_idx));
    idx_contain_pk = []; 
    for j = 1:length(run_inds)
        if any(run_inds{j} == p_locs(valid_peak_idx))
            idx_contain_pk = j;
            break;
        end
    end
    exc_period = run_len(idx_contain_pk)*xbin;
end

% Save results
hist_features = struct();
hist_features.ref_period = ref_period;
hist_features.exc_period = exc_period;
hist_features.p_time = p_time;
hist_features.p_height = p_height;

% Save infraslow multiplier only when 
if xlag(end) > 1500
    is_left = 40/xbin;
    is_right = 70/xbin;
    is_length = (is_right-is_left)*xbin;
    AUC_is = sum((yhat(is_left:is_right))*xbin)/is_length;
    hist_features.AUC_is = AUC_is;
end

end
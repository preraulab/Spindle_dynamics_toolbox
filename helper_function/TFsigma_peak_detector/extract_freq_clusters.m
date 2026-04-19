function [ freq_TFpeaks ] = extract_freq_clusters(freq_TFpeaks, sel_freqs)
%EXTRACT_FREQ_CLUSTERS  Identify TFpeak indices that fall within the extent of each frequency cluster
%
%   Usage:
%       freq_TFpeaks = extract_freq_clusters(freq_TFpeaks, sel_freqs)
%
%   Inputs:
%       freq_TFpeaks : table - frequency cluster table with peak_lower_freq, peak_upper_freq, and boundary_from_lastpeak columns -- required
%       sel_freqs    : 1xN double - central frequencies of candidate TFpeaks -- required
%
%   Outputs:
%       freq_TFpeaks : table - input augmented with a TFpeak_idx column (cell of logical vectors)
%
%   Notes:
%       Low/high extent per cluster is widened by half the cluster width,
%       and clamped by the trough between adjacent clusters.
%
%   See also: extract_maxfreq_peaks, extract_density_curve
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

for ii = 1:size(freq_TFpeaks,1)
    low_bound = freq_TFpeaks.peak_lower_freq(ii) - (freq_TFpeaks.peak_upper_freq(ii) - freq_TFpeaks.peak_lower_freq(ii))/2;
    low_bound = max(freq_TFpeaks.boundary_from_lastpeak(ii), low_bound);
    
    high_bound = freq_TFpeaks.peak_upper_freq(ii) + (freq_TFpeaks.peak_upper_freq(ii) - freq_TFpeaks.peak_lower_freq(ii))/2;
    if ii ~= size(freq_TFpeaks,1)
        high_bound = min(freq_TFpeaks.boundary_from_lastpeak(ii+1), high_bound);
    end
    
    % identify the indices 
    peak_indices = sel_freqs >= low_bound & sel_freqs <= high_bound;
    freq_TFpeaks.TFpeak_idx(ii) = {peak_indices};
end


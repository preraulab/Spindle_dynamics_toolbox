function [ max_curve, bin_centers, hist_olN2, hist_olN3 ] = extract_density_curve(sel_freqs, sel_stages, bin_width, bin_step, bin_range, N2_minutes, N3_minutes, ignore_N3_threshold, plot_on)
%EXTRACT_DENSITY_CURVE  Extract per-minute frequency density curves for N2/N3 TFpeaks and their max
%
%   Usage:
%       [max_curve, bin_centers, hist_olN2, hist_olN3] = extract_density_curve( ...
%           sel_freqs, sel_stages, bin_width, bin_step, bin_range, ...
%           N2_minutes, N3_minutes, ignore_N3_threshold, plot_on)
%
%   Inputs:
%       sel_freqs           : 1xN double - central frequencies (Hz) of selected TFpeaks -- required
%       sel_stages          : 1xN categorical - sleep stage of each TFpeak -- required
%       bin_width           : double - frequency bin width in Hz -- required
%       bin_step            : double - frequency bin step in Hz -- required
%       bin_range           : 1x2 double - [min, max] frequency range for binning -- required
%       N2_minutes          : double - total N2 duration in minutes (for density normalization) -- required
%       N3_minutes          : double - total N3 duration in minutes (for density normalization) -- required
%       ignore_N3_threshold : double - ignore N3 curve if N3_minutes < this threshold -- required
%       plot_on             : double - 0 = no plot, 1 = new figure, otherwise handle to target axes -- required
%
%   Outputs:
%       max_curve   : 1xB double - per-bin max of the N2 and N3 density curves (events/min)
%       bin_centers : 1xB double - centers of frequency bins (Hz)
%       hist_olN2   : 1xB double - N2 binned density curve (events/min)
%       hist_olN3   : 1xB double - N3 binned density curve (events/min)
%
%   Notes:
%       Heuristic of max(stage2, stage3) is used to build a combined
%       density curve for mid-point frequencies of TFpeaks.
%
%   See also: extract_maxfreq_peaks, extract_freq_clusters
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

% set up bin starts and ends
hb = bin_width/2;
bin_centers = bin_range(1)+hb:bin_step:bin_range(2)-hb;
bin_starts = bin_centers - hb;
bin_ends = bin_centers + hb;

% stage 2 binned cumulative density
data = sel_freqs(sel_stages=='Stage2');
hist_ol = zeros(1, length(bin_centers));
for jj = 1:length(bin_centers)
    hist_ol(jj) = sum(data>bin_starts(jj) & data<=bin_ends(jj));
end
hist_olN2 = hist_ol / N2_minutes;

% stage 3 binned cumulative density
data = sel_freqs(sel_stages=='Stage3');
hist_ol = zeros(1, length(bin_centers));
if N3_minutes > ignore_N3_threshold % ignore stage 3 if duration of N3 is shorter than [ignore_N3_threshold] minutes
    for jj = 1:length(bin_centers)
        hist_ol(jj) = sum(data>bin_starts(jj) & data<=bin_ends(jj));
    end
end
hist_olN3 = hist_ol / N3_minutes;

if plot_on ~= 0
    if plot_on == 1
        figure
    else
        axes(plot_on)
    end
    hold on
    plot(bin_centers, hist_olN2, 'Linewidth', 2)
    plot(bin_centers, hist_olN3, 'Linewidth', 2)
    legend(['Stage2=', num2str(N2_minutes), 'min'], ['Stage3=', num2str(N3_minutes), 'min'], 'Location', 'best')
    xlabel('Frequency (Hz)')
    ylabel('Binned sum density (events / min)')
    title('Binned sum density')
    set(gca,'FontSize',16)
    xlim([min(bin_centers), max(bin_centers)])
    xticks([min(bin_centers):1:max(bin_centers)])
end

% extract the max out of the two curves
max_curve = max([hist_olN2; hist_olN3]);

end


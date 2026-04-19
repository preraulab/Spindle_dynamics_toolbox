%EXAMPLE_SCRIPT  Run the full spindle-dynamics analysis pipeline on example data and generate summary figures
%
%   Usage:
%       Example_Script
%
%   Inputs:
%       none (loads example_data/example_data1.mat from the repository)
%
%   Outputs:
%       none (side effects only - produces analysis figures)
%
%   Notes:
%       Walks through model specification (specify_mdl), preprocessing
%       (preprocessToDesignMatrix), Poisson GLM fitting (glmfit), and the
%       full set of summary plots (history curve, stage preferred phase,
%       SO-power preferred phase, KS goodness-of-fit, deviance explained).
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: quick_start, specify_mdl, preprocessToDesignMatrix, plot_hist_curve, plot_stage_prefphase, plot_sop_prefphase
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% PATH SETTINGS

% Add help functions to path
addpath(genpath('./helper_function'))

%% DATA LOADING
% Location of example data
data_path = 'example_data/';

% Load example EEG data 1
load([data_path,'example_data1.mat']);

%% PRE SETTINGS AND MODEL SPECIFICATIONS
binsize = 0.1;                  % point process bin size in sec (Fixed)
hard_cutoffs = [12, 16];        % choose fast spindles in freq range:12-16 Hz

% 1) Choose model factors by a binary vector (1x4,double)
%    in the fixed order: SOphase, stage, SOpower,history
%    This vector determines which factors are included (1) or not (0)

BinarySelect = [1,1,0,1];  % Select SOphase,stage,history as factors

% 2) Choose interaction terms in the form of a cell {'XxX:YYy'}   
%    It is case,order-insensitive,and accepts separators such as ':', '&', and '-'

InteractSelect = {'stage:SOphase'};  % Add stage-SOphase interaction

% 3) Choose history order
%   'long': Long term history (up to 90 secs, this can show infraslow activity)
%   'short': Short term history (up to 15 secs, this option runs fast)

hist_choice = 'long';
switch lower(hist_choice)
    case {'short'}
          hist_ord = 15/binsize;   % 15 sec history
          control_pt = [0:15:90 120 hist_ord]; % cardinal spline control point location
    case {'long'}
          hist_ord = 90/binsize;   % 90 sec history
          control_pt = [0:15:90 120 150:100:750 hist_ord]; 
end

% 4) Specify the model
[ModelSpec] = specify_mdl(BinarySelect,InteractSelect,'hist_choice',hist_choice,'control_pt',control_pt);

%% DATA PREPROCESS TO DESIGN MATRIX
[X, BinData,res_table] = preprocessToDesignMatrix(EEG, Fs, stage_val, stage_time,ModelSpec);

%% MODEL FITTING
warning('off', 'all');
[b, dev, stats] = glmfit(X,BinData.y,'poisson'); 

%% PREPARE FOR FIGURE
% Compute multitaper spectrogram
dsfreqs = 0.1;
nfft = 2^(nextpow2(Fs/dsfreqs));
[mt, stimes, sfreqs] = multitaper_spectrogram_mex(EEG, Fs,[8,19],[2,3],[1,0.05],nfft,'constant',[],false,false);

% Extract fast spindle raw time and freq
fast_inds = (res_table.peak_freqs{:} >= hard_cutoffs(1)) & (res_table.peak_freqs{:}< hard_cutoffs(2));
spindletime_all = res_table.peak_ctimes{:};
freq_all  = res_table.peak_freqs{:};
fast_pt = spindletime_all(fast_inds);
fast_pf = freq_all(fast_inds);

% Extract raw SO power and phase
phase_raw = cell2mat(res_table.SOphase);      % SO phase raw
sop_raw = cell2mat(res_table.SOpow);          % SO power raw

% Extract data from BinData
xtime = BinData.real_time + binsize/2;        % bin data time stamps 
sta = BinData.sta;                            % stage binary matrix [N1,N2,N3,REM,Wake]  
stagetrain = 3*sta(:,1)+ 2*sta(:,2)+ sta(:,3)+4*sta(:,4)+5*sta(:,5); % stagetrain in bin
y = BinData.y;                                % spindle binary train
sop = BinData.sop;                            % SO power in bin  
phase = BinData.phase;                        % SO phase in bin  

%% GENERATE THE OVERVIEW FIGURE
time_range = [xtime(1) xtime(end)];  % Set the time range of the visualization

% Figure design
figure; 
ax = figdesign(7, 3, 'orient', 'portrait' , 'margins', [0.08 0.11 0.08 0.08 0.08 0.08], ...
           'merge', {[1 2 4 5], [7 8], [10 11], [13 14], [16 17], [19 20], [3 6 9 12], [15 18 21]}, ...
           'Position', [0.01 0.01 1 1]);
linkaxes(ax([1 2 3 5 6 7]),'x')

% ---------------------------Left panel ---------------
% 1) Spectrogram
axes(ax(1))
hold on;
imagesc(stimes, sfreqs, nanpow2db(mt));
axis xy;
ylabel('Frequency (Hz)');
climscale;
c = colorbar_noresize;
colormap(gca,rainbow4);
ylabel(c,'Power (dB)');
axis tight
plot(fast_pt,fast_pf,'o','markersize', 6,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1.5);
title('Sleep EEG Spectrogram')
scaleline(ax(1),30,'30 secs');

% 2) Spindle train
axes(ax(2))
stem(fast_pt,ones(size(fast_pt)),'Color','k', 'Marker', 'none','linewidth',1.5);
set(ax(2),'YColor','none')
ylim(ax(2),[0 1.3])
title('Sleep Spindle Train')

% 3) Stage train
axes(ax(3))
hold on;
hypnoplot(stage_time, stage_val);
plot(xtime(y==1),stagetrain(y==1),'o','markersize', 6,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1.5);
title('Sleep Stage')

% 4) SO power
axes(ax(5))
hold on;
plot(t,sop_raw,'LineWidth',1.5,'Color','k');
plot(xtime(y==1),sop(y==1),'o','markersize', 6,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1.5);
ylabel('Power')
title('Continuous Sleep Depth (SO Power)');

% 5) SO phase
axes(ax(6))
hold on;
plot(t,phase_raw,'k','LineWidth',1.5);
phase  = wrapToPi(phase);
plot(xtime(y==1),phase(y==1),'o','markersize', 6,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1.5);
set(ax(6),'YTick',[-pi 0 pi],'YTickLabels',{'-\pi', '0', '\pi'},'ylim',[-pi pi])
yline(0,'k--')
ylabel('Phase')
title('Cortical Up/Down State (SO Phase)');

% 6) CIF
axes(ax(7))
[CIF_sta,~,~] = glmval(b,X,'log',stats); % Compute instantaneous spindle rate
hold on;
plot(xtime,CIF_sta*60/binsize,'color','k','LineWidth',2);
plot(xtime(y==1),CIF_sta(y==1)*60/binsize,'o','markersize', 6,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',2);
title(sprintf('Density = f (%s)', ModelSpec.FactorSelect));

xlabel('Time (s)',FontSize=12);
ylabel({'Spindle Density';'(events/min)'},FontSize=12);
pos = get(ax(7), 'Position'); 
pos(4) = pos(4) + 0.05; 
pos(2) = pos(2) - 0.05; 
set(ax(7), 'Position', pos); 
xlim(time_range)
scaleline(ax(7),5,'5 secs');

scrollzoompan(ax(7),'x');   
% To get the same overview figure in the tutorial:
% In Zoom box, input 90, hit enter
% In Pan box, input 30145, hit enter 


%--------------------- Right panel -----------------------
% Plot Polarhistogram
axes(ax(4))
phasehistogram(phase(y==1),1,'NumBins',50); 
rlim([0 .4])
rticks(0:0.1:0.4)
title('Spindle Phase Histogram',FontSize=12)

% Plot history Curve
axes(ax(8))
plot_hist_curve(stats,ModelSpec,BinData);

% Fig settings
for i = [1 2 3 5 6 7]
    title(ax(i),ax(i).Title.String,'FontSize',12,'FontWeight','bold'); 
    set(ax(i),'FontSize',12,'box','off');
end

%% GENERATE THE HISTORY MODULATION FIGURE
figure;
[xlag,yhat,yu,yl,hist_features] = plot_hist_curve(stats,ModelSpec,BinData);

%% GENERATE THE STAGE VS SOP PHASE SHIFT FIGURE
% Fit model with SOP
[ModelSpec_sop] = specify_mdl([1 0 1 1],{'SOphase:SOpower'},'hist_choice',hist_choice,'control_pt',control_pt);
X_sop = build_design_mt(ModelSpec_sop,BinData);
[b_sop, dev_sop, stats_sop] = glmfit(X_sop,y,'poisson'); 
[CIF_sop,~,~] = glmval(b_sop,X_sop,'log',stats_sop);

% Figure
figure; 
ax = figdesign(4,6,'margins',[.08 .13 .02 .08 .05 .05],'merge',{[2 3 8 9 14 15],[4 10 16],[5 6 11 12 17 18],[20 21],[23 24]});
set(gcf, 'units','inches','Position',  [0, 0, 12, 5]);
linkaxes(ax([5 6]),'y')
delete(ax(7))
delete(ax(9))

%----------- Stage PP polarhistogram
% N1
axes(ax(1))
[~, ~, h_phist, h_pax, ~] = phasehistogram(phase(sta(:,1)==1&y==1),1,'NumBins',50); 
rlim([0 .4])
rticks(0:0.1:0.4)
title('N1',FontSize=12)

% N2
axes(ax(2))
phasehistogram(phase(sta(:,2)==1&y==1),1,'NumBins',50); 
rlim([0 .4])
rticks(0:0.1:0.4)
title('N2',FontSize=12)

% N3
axes(ax(3))
phasehistogram(phase(sta(:,3)==1&y==1),1,'NumBins',50); 
rlim([0 .4])
rticks(0:0.1:0.4)
title('N3',FontSize=12)

%----------- Stage PP shift 
axes(ax(4))
[pp,pp_CI] = plot_stage_prefphase(b,stats,ModelSpec);
title('Model With Stage')

%----------- Stage percentage
axes(ax(5))
[sta_prop,sop_xvalue] = compute_stage_prob(stagetrain, sop, linspace(min(sop),max(sop),100)'); 
pn1 = plot(sta_prop(:,3)*100,sop_xvalue,'b','linewidth',2.5);hold on;axis ij;
pn2 = plot(sta_prop(:,2)*100,sop_xvalue,'r','linewidth',2.5);
pn3 = plot(sta_prop(:,1)*100,sop_xvalue,'g','linewidth',2.5);
xlabel('% Time in Stage','fontsize',10);
ylabel('SO Power','fontsize',12);
lg = legend([pn1 pn2 pn3],{'N1','N2','N3'},'fontsize',10,'Location','northeast','orientation','vertical','Box','off','Color','none');
lg.ItemTokenSize = [7, 7];
ax1 = gca; 
ax1.FontSize = 12;
positionAdj = get(ax(5), 'Position'); % Get and adjust the position
positionAdj(1) = positionAdj(1) + 0.02; % Move right by 0.05 (adjust as needed)
set(ax(5), 'Position', positionAdj);% Set the new position

%----------- SOP PP shift curve
axes(ax(6))
[phi0,sop0,sop_pp_mat] = plot_sop_prefphase(b_sop,stats_sop,ModelSpec_sop);
ylim([3 25])
set(gca,'YColor','none')

% Set ratemaps on the same scale
clim_sta = clim(ax(4)); 
clim_sop = clim(ax(6)); 
comClim = [min(clim_sta(1), clim_sop(1)), max(clim_sta(2), clim_sop(2))];
clim(ax(4),comClim)
clim(ax(6),comClim)

% ----------- Add cosine reference plots
axes(ax(8))
plotcos;
axes(ax(10))
plotcos;

%% GENERATE KS PLOTS (MODEL EVALUATION)
% Stage only model
[ModelSpec1_s] = specify_mdl([0 1 0 0],[],'hist_choice',[],'control_pt',[]);
X1_s = build_design_mt(ModelSpec1_s,BinData);
[b1_s, dev1_s, stats1_s] = glmfit(X1_s,y,'poisson'); 
[CIF1_s,~,~] = glmval(b1_s,X1_s,'log',stats1_s);

% Phase only Model
[ModelSpec1_p] = specify_mdl([1 0 0 0],[],'hist_choice',[],'control_pt',[]);
X1_p = build_design_mt(ModelSpec1_p,BinData);
[b1_p, dev1_p, stats1_p] = glmfit(X1_p,y,'poisson'); 
[CIF1_p,~,~] = glmval(b1_p,X1_p,'log',stats1_p);

% History only model
[ModelSpec1_h] = specify_mdl([0 0 0 1],[]);
X1_h = build_design_mt(ModelSpec1_h,BinData);
[b1_h, dev1_h, stats1_h] = glmfit(X1_h,y,'poisson'); 
[CIF1_h,~,~] = glmval(b1_h,X1_h,'log',stats1_h);

%--------- KS Figure
figure; 
ax = figdesign(1,4,'margins',[.05 .05 .05 .05 .07 .05]);
set(gcf, 'units','inches','Position',  [0, 0, 14, 4]);

axes(ax(1))
[ks1_s,ksT1_s] =  KSplot(CIF1_s,y,1);  % 0 means pass the KS test
title('Model with only stage')

axes(ax(2))
[ks1_p,ksT1_p] =  KSplot(CIF1_p,y,1);
title('Model with only phase')

axes(ax(3))
[ks1_h,ksT1_h] =  KSplot(CIF1_h,y,1); 
title('Model with only history')

axes(ax(4))
[ks_s,ksT_s] =  KSplot(CIF_sta,y,1); 
title({'Model with stage, phase, history','and stage-phase interaction'})

for i = 1:4
    axis(ax(i), 'square'); 
end

%% OUTPUT DEVIANCE EXPLAINED FOR EACH FACTOR
[dev_exp_sta,dev_exp_sop] = compute_dev_exp(BinData);

% Display the table
factors = {'Stage'; 'Phase'; 'History'};
Table1 = table(factors, 100*dev_exp_sta', 'VariableNames', {'Factor', 'Deviance Explained %'});
disp(Table1);


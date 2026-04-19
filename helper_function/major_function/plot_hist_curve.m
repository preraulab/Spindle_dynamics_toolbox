function [xlag,yhat,yu,yl,hist_features] = plot_hist_curve(stats,ModelSpec,BinData)
%PLOT_HIST_CURVE  Compute history modulation curves, derived features, and plot the history curve
%
%   Usage:
%       [xlag, yhat, yu, yl, hist_features] = plot_hist_curve(stats, ModelSpec, BinData)
%
%   Inputs:
%       stats     : struct - GLM results returned by glmfit -- required
%       ModelSpec : struct - model specification (binsize, hist_ord, control_pt, ...) -- required
%       BinData   : struct - binned data (sp, isis, ...) -- required
%
%   Outputs:
%       xlag          : nx1 double - history time lag in seconds
%       yhat          : nxk double - history modulation (rate multiplier) values
%       yu            : nxk double - upper 95 percent confidence bound
%       yl            : nxk double - lower 95 percent confidence bound
%       hist_features : struct - derived features of the history curve
%
%   Notes:
%       k = 1 for a single history curve; k = 2 when stage-dependent N2/N3
%       curves are requested (col 1 = N2, col 2 = N3). n = hist_ord / sp_resol.
%       hist_features fields:
%           ref_period : refractory period (s)
%           exc_period : excited period (s)
%           p_time     : peak time (s)
%           p_height   : peak height
%           AUC_is     : area under the infraslow period (40-70 s), only when
%                        'long' history is specified
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: compute_hist_features, FinerModCardinalSpline, glmval, shadebounds
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

if ModelSpec.BinarySelect(4) == 1
%% Prepare for the figure
b = stats.beta;              % fitted parameters 
sp_resol = 0.1;              % Evaluate spline in a finer resolution
[sp_finer] = FinerModCardinalSpline(ModelSpec.hist_ord,ModelSpec.control_pt,.5,sp_resol);

% time lag in sec
xlag = (ModelSpec.binsize*sp_resol:ModelSpec.binsize*sp_resol:ModelSpec.hist_ord*ModelSpec.binsize)';

    %% Evaluate history modulation curve
    %--- If single history curve is specified
    if length(ModelSpec.InteractSelect)<=1
        sizediff = length(b)-1 - size(BinData.sp,2); % num of non-history columns in design matrix
        [yhat0,yl0,yu0] = glmval(b,[zeros(ModelSpec.hist_ord/sp_resol,sizediff) sp_finer],'log',stats);
        yl = (yhat0 - yl0)/exp(b(1));
        yu = (yhat0 + yu0)/exp(b(1));
        yhat = yhat0/exp(b(1));
        [hist_features] = compute_hist_features(xlag,yhat,yl,yu);

        % figure
        shadebounds(xlag,yhat,yu,yl,'k',[.5,.5,.5],[.9,.9,.9],.4);
        stem(BinData.isis, -0.2*ones(length(BinData.isis),1),'Color','k', 'Marker', 'none','linewidth',.05);
        yline(1,'k--')
        xlim([0 ModelSpec.hist_ord*ModelSpec.binsize])
        ylim([-0.2 2])
        xlabel('Time Lag (sec)')
        ylabel('Rate Multiplier')
        title('History Modulation');
        set(gca,'XTick',0:5:ModelSpec.hist_ord*ModelSpec.binsize,'fontsize',12)
    
    %--- If stage-dependent history is specified
    elseif length(ModelSpec.InteractSelect)==2
        sizediff = length(b)-1 - 3*size(BinData.sp,2); % num of non-history columns in design matrix
        % Compute N2 history curve
        [y2,yl2,yh2] = glmval(b,[zeros(ModelSpec.hist_ord/sp_resol,sizediff) sp_finer sp_finer zeros(size(sp_finer))],'log',stats);
        yl2 = (y2 - yl2)/exp(b(1));
        yu2 = (y2 + yh2)/exp(b(1));
        y2 = y2/exp(b(1));
        [hist_features] = compute_hist_features(xlag,y2,yl2,yu2);
        
        % Compute N3 history curve
        [y3,yl3,yh3] = glmval(b,[zeros(ModelSpec.hist_ord/sp_resol,sizediff) sp_finer zeros(size(sp_finer)) sp_finer],'log',stats);
        yl3 = (y3 - yl3)/exp(b(1));
        yu3 = (y3 + yh3)/exp(b(1));
        y3 = y3/exp(b(1));
    
        % Save results
        yl = [yl2 yl3];
        yu = [yu2 yu3];
        yhat = [y2 y3];
    
        % figure
        hold on;
        s1 = shadebounds(xlag,y2,yu2,yl2,'b','b','none',.2);
        s2 = shadebounds(xlag,y3,yu3,yl3,'m','m','none',.2);
        stem(BinData.isis, -0.2*ones(length(BinData.isis),1),'Color','k', 'Marker', 'none','linewidth',.05);
        yline(1,'k--')
        xlim([0 ModelSpec.hist_ord*ModelSpec.binsize])
        ylim([-0.2 2])
        xlabel('Time Lag (sec)')
        ylabel('Rate Multiplier')
        title('History Modulation');
        legend([s1,s2],{'N2 History Curve','N3 History Curve'})
        set(gca,'XTick',0:5:ModelSpec.hist_ord*ModelSpec.binsize,'fontsize',12)
    end

else
    error('History component is not selected by the user')
end


end
function [phi0,sop0,sop_pp_mat] = plot_sop_prefphase(b,stats,ModelSpec)
%PLOT_SOP_PREFPHASE  Plot preferred SO phase as a function of continuous SO power
%
%   Usage:
%       [phi0, sop0, sop_pp_mat] = plot_sop_prefphase(b, stats, ModelSpec)
%
%   Inputs:
%       b         : 1xM double - fitted model parameters from glmfit -- required
%       stats     : struct - glmfit stats (covb, ...) -- required
%       ModelSpec : struct - model specification, must contain SOpower:SOphase interaction -- required
%
%   Outputs:
%       phi0       : 500x1 double - phase samples, 500 evenly spaced points in [-pi, pi]
%       sop0       : 510x1 double - SOP samples, 510 evenly spaced points
%       sop_pp_mat : 510x500 double - spindle rate heatmap as a function of (SOP, phase)
%
%   Notes:
%       Requires SOpower:SOphase to be in ModelSpec.InteractSelect. The 95
%       percent confidence band for the preferred phase is estimated by
%       Monte Carlo sampling (M = 5000) from the coefficient covariance.
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: plot_stage_prefphase, specify_mdl, ismember_interaction
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Evaluate rate heatmap ~ phase,sop from the model
% SOpower:SOphase Must be specified to generate this
if any(ismember_interaction(ModelSpec.InteractSelect,{'SOpower:SOphase'}))

    phi0 = linspace(-pi,pi,500)';
    sop0 = linspace(-1,6,510)';
    scale_factor = 10;
    
    % Initialization
    sop_pp_mat = zeros(length(sop0),length(phi0));
    if ModelSpec.BinarySelect(4) == 0
        sizeh = [];
    else
        sizeh = zeros(length(phi0),length(ModelSpec.control_pt));
    end

    for j = 1:length(sop0)
           [y_sop, ~, ~] = glmval(b,[zeros(size(phi0)) ones(size(phi0))*sop0(j) ones(size(phi0)).*(sop0(j)^2)/scale_factor cos(phi0) sin(phi0)  ...
                                    cos(phi0).*sop0(j) cos(phi0).*(sop0(j)^2)/scale_factor ...
                                    sin(phi0).*sop0(j) sin(phi0).*(sop0(j)^2)/scale_factor sizeh],'log',stats);
           sop_pp_mat(j,:) = y_sop;
    end    
    
    %% Pref phase ~ sop
    pp_sop = zeros(length(sop0),1);
    
    for i = 1:length(sop0)
         pp_sop(i) = atan2(b(6)+b(9)*sop0(i) + b(10)*sop0(i)^2/scale_factor, ...
                           b(5)+b(7)*sop0(i) + b(8)*sop0(i)^2/scale_factor); 
    end
    
    %% Confidence interval for SOP PP shift curve
    M = 5000;          % sampling size, 10000 is used in the paper
    bMC = mvnrnd(b,stats.covb,M)';
    pp_sop_lohi = zeros(length(sop0),2);
    
    for j = 1:length(sop0)
    
         yMC = glmval(bMC,[zeros(size(phi0)) ones(size(phi0))*sop0(j) ones(size(phi0)).*(sop0(j)^2)/scale_factor cos(phi0) sin(phi0)  ...
                           cos(phi0).*sop0(j) cos(phi0).*(sop0(j)^2)/scale_factor ...
                           sin(phi0).*sop0(j) sin(phi0).*(sop0(j)^2)/scale_factor sizeh],'log',stats);
         mx = zeros(M,1);                                                             
         for k = 1:M
             [~,mx1]  =  max(yMC(:,k));
             mx(k) = phi0(mx1);
         end
         
         pp_sop_lohi(j,:) = quantile(mod(mx-pp_sop(j) + pi,2*pi)-pi, [0.025, 0.975]);
    
    end
    
    % Save CI
    pp_CI = [pp_sop_lohi(:,1) + pp_sop pp_sop_lohi(:,2) + pp_sop];
    
    %% Prepare for figure
    % Choose not to plot CI if it crosses [-pi pi]
    pp_CI_plot = pp_CI;
    
    % for i = 1:size(pp_CI,1)
    %     if pp_CI(i,1)<-pi||pp_CI(i,2)>pi
    %         pp_CI_plot(i,1) = -pi-0.5;
    %         pp_CI_plot(i,2) =  pi+0.5;
    %     end
    % end
    
    %% Figure
    sop0 = sop0*5;
    imagesc(phi0,sop0,sop_pp_mat*60/ModelSpec.binsize); % spindle rate heatmap, transform to events/min
    hold on;
    colormap(gca,viridis);
    plt_pp = plot(pp_sop,sop0,'r',LineWidth = 3);
    plt_ci = plot(pp_CI_plot(:,1),sop0,'r--');
             plot(pp_CI_plot(:,2),sop0,'r--');
    xlim([-pi pi])
    
    set(gca,'YTick',[],'XTick',[])
    a2 = colorbar_noresize;
    a2.Location = 'Eastoutside';
    ylabel(a2,'Rate (events/min)','FontSize',12);
    clim([0 max(sop_pp_mat'*60/ModelSpec.binsize,[],'all')])
    title('Model with SOP','Fontsize',12);
    legend([plt_pp plt_ci],{'Preferred Phase','Confidence Interval'},'fontsize',10,'Location','northeast');
    ax1 = gca; 
    ax1.FontSize = 12;
    hold off;
    
else 
    error('This can not be evaluated as SOpower-SOphase interaction is not included in the model')

end
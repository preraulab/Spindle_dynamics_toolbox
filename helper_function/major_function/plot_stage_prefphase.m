function [pp,pp_CI] = plot_stage_prefphase(b,stats,ModelSpec)
%PLOT_STAGE_PREFPHASE  Plot preferred SO phase in N1, N2, and N3 sleep stages
%
%   Usage:
%       [pp, pp_CI] = plot_stage_prefphase(b, stats, ModelSpec)
%
%   Inputs:
%       b         : 1xM double - fitted model parameters from glmfit -- required
%       stats     : struct - glmfit stats (covb, ...) -- required
%       ModelSpec : struct - model specification, must contain stage:SOphase interaction -- required
%
%   Outputs:
%       pp    : 3x1 double - preferred phase (rad) for N1, N2, N3
%       pp_CI : 3x2 double - 95 percent CI for each preferred phase (col 1: lower, col 2: upper)
%
%   Notes:
%       Requires stage:SOphase to be in ModelSpec.InteractSelect. The CI is
%       estimated by Monte Carlo sampling (M = 10000) from the coefficient
%       covariance. Accompanies Chen et al., PNAS 2025.
%
%   See also: plot_sop_prefphase, specify_mdl, ismember_interaction
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% stage:SOphase must be specified to generate this
if any(ismember_interaction(ModelSpec.InteractSelect,{'stage:SOphase'}))
    %% Compute preferred phase, magnitude from the model
    
    pp_n1 = atan2(b(3)+b(12),b(2)+b(8)); % pref phase for N1
    pp_n2 = atan2(b(3),b(2));            % pref phase for N2
    pp_n3 = atan2(b(3)+b(13),b(2)+b(9)); % pref phase for N3
    pp = [pp_n1; pp_n2; pp_n3];            % Save outputs

    % mag_n1 = sqrt((b(3)+b(12))^2 + (b(2)+b(8))^2); % coupling magnitude for N1
    % mag_n2 = sqrt(b(3)^2 + b(2)^2);                % coupling magnitude for N2
    % mag_n3 = sqrt((b(3)+b(13))^2 + (b(2)+b(9))^2); % coupling magnitude for N3

    %% Evaluate rate ~ phase from the model
    % Initialization
    phi0 = linspace(-pi,pi,500)';
    if ModelSpec.BinarySelect(4) == 0
        sizeh = [];
    else
        if all(ismember_interaction(ModelSpec.InteractSelect,{'stage:SOphase'}))
            sizeh = zeros(length(phi0),length(ModelSpec.control_pt));
        elseif all(ismember_interaction(ModelSpec.InteractSelect,{'stage:SOphase','stage:history'}))
            sizeh = zeros(length(phi0),3*length(ModelSpec.control_pt));
        end
    end

    % N1
    [y_n1, y_n1_hi,y_n1_lo] = glmval(b,[cos(phi0) sin(phi0) ones(length(phi0),1) zeros(length(phi0),3) ...
                                        cos(phi0).*ones(length(phi0),1) cos(phi0).*zeros(length(phi0),3) ...
                                        sin(phi0).*ones(length(phi0),1) sin(phi0).*zeros(length(phi0),3) sizeh],'log',stats);
    % N2
    [y_n2, y_n2_hi,y_n2_lo] = glmval(b,[cos(phi0) sin(phi0) zeros(length(phi0),12) sizeh],'log',stats);
    
    % N3
    [y_n3, y_n3_hi,y_n3_lo] = glmval(b,[cos(phi0) sin(phi0) zeros(length(phi0),1) ones(length(phi0),1) zeros(length(phi0),2)...
                                        cos(phi0).*zeros(length(phi0),1) cos(phi0).*ones(length(phi0),1) cos(phi0).*zeros(length(phi0),2)...
                                        sin(phi0).*zeros(length(phi0),1) sin(phi0).*ones(length(phi0),1) sin(phi0).*zeros(length(phi0),2) sizeh],'log',stats);
    %% Compute confidence interval of the preferred phase
    M = 10000;          % sampling size
    bMC = mvnrnd(b,stats.covb,M)';
    
    % N1
    yMC_n1 = glmval(bMC,[cos(phi0) sin(phi0) ones(length(phi0),1) zeros(length(phi0),3)...
                      cos(phi0).*ones(length(phi0),1) cos(phi0).*zeros(length(phi0),3)...
                      sin(phi0).*ones(length(phi0),1) sin(phi0).*zeros(length(phi0),3) sizeh],'log',stats);
    mx = zeros(M,1);
    for k=1:M
      [~,mx1]  =  max(yMC_n1(:,k));
      mx(k) = phi0(mx1);
    end
    pp_n1_lohi = quantile(mod(mx-pp_n1 + pi,2*pi)-pi, [0.025, 0.975]);
    
    % N2
    yMC_n2 = glmval(bMC,[cos(phi0) sin(phi0) zeros(length(phi0),12) sizeh],'log',stats);
    mx = zeros(M,1);
    for k=1:M
      [~,mx2]  =  max(yMC_n2(:,k));
      mx(k) = phi0(mx2);
    end
    pp_n2_lohi = quantile(mod(mx-pp_n2+pi,2*pi)-pi, [0.025, 0.975]);
    
    %N3
    yMC_n3 = glmval(bMC,[cos(phi0) sin(phi0) zeros(length(phi0),1) ones(length(phi0),1) zeros(length(phi0),2)...
                         cos(phi0).*zeros(length(phi0),1) cos(phi0).*ones(length(phi0),1) cos(phi0).*zeros(length(phi0),2)...
                         sin(phi0).*zeros(length(phi0),1) sin(phi0).*ones(length(phi0),1) sin(phi0).*zeros(length(phi0),2) sizeh],'log',stats);
    mx = zeros(M,1);
    for k=1:M
      [~,mx3]  =  max(yMC_n3(:,k));
      mx(k) = phi0(mx3);
    end
    pp_n3_lohi = quantile(mod(mx-pp_n3+pi,2*pi)-pi, [0.025, 0.975]);
    
    % Save CI for N1 pref phase, N2 pref phase, and N3 pref phase
    pp_CI = [pp_n1_lohi + pp_n1; pp_n2_lohi + pp_n2; pp_n3_lohi + pp_n3];
    
    %% Prepare for figure
    % Choose not to plot CI if it crosses [-pi pi]
    pp_CI_plot = pp_CI;
    for i = 1:size(pp_CI,1)
        if pp_CI(i,1)<-pi||pp_CI(i,2)>pi
            pp_CI_plot(i,1) = -pi-0.5;
            pp_CI_plot(i,2) =  pi+0.5;
        end
    end
    
    % Generate a stage pp matrix: Rate ~ phase heatmap
    stage_pp_mat = zeros(length(phi0),600);
    n1_cut = [1 200];
    n2_cut = [200 400];
    n3_cut = [400 600];
    stage_pp_mat(:,1:200)   = repmat(y_n1,1,200);
    stage_pp_mat(:,201:400) = repmat(y_n2,1,200);
    stage_pp_mat(:,401:600) = repmat(y_n3,1,200);

    
    %% Figure
    
    imagesc(phi0,1:600,stage_pp_mat'*60/ModelSpec.binsize); hold on;
    colormap(gca,viridis);
    pp_line = line([pp_n1 pp_n1],n1_cut,'color','r','linestyle','-','linewidth',3);
              line([pp_n2 pp_n2],n2_cut,'color','r','linestyle','-','linewidth',3);
              line([pp_n3 pp_n3],n3_cut,'color','r','linestyle','-','linewidth',3);
    %lower bound
    pp_line_ci = line([pp_CI_plot(1,1) pp_CI_plot(1,1)],n1_cut,'color','r','linestyle','--');
                 line([pp_CI_plot(2,1) pp_CI_plot(2,1)],n2_cut,'color','r','linestyle','--');
                 line([pp_CI_plot(3,1) pp_CI_plot(3,1)],n3_cut,'color','r','linestyle','--');
    %upper bound
    line([pp_CI_plot(1,2) pp_CI_plot(1,2)],n1_cut,'color','r','linestyle','--');
    line([pp_CI_plot(2,2) pp_CI_plot(2,2)],n2_cut,'color','r','linestyle','--');
    line([pp_CI_plot(3,2) pp_CI_plot(3,2)],n3_cut,'color','r','linestyle','--');
    
    legend([pp_line pp_line_ci],{'Preferred Phase','Confidence Interval'},'fontsize',10,'Location','northwest');
    title('Preferred Phase From The Model','Fontsize',12);
    xlim([-pi pi])
    yticks([100 300 500]);
    yticklabels({'N1','N2','N3'});
    set(gca,'XTick',[])
    ax1 = gca; 
    ax1.FontSize = 12;
    hold off
  
else 
    error('This can not be evaluated as stage-SOphase interaction is not included in the model')

end
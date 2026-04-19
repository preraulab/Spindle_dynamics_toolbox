function [] = plotcos()
%PLOTCOS  Draw a cosine reference plot annotated with SO up/down states
%
%   Usage:
%       plotcos()
%
%   Inputs:
%       none
%
%   Outputs:
%       none (side effects only - draws into the current axes)
%
%   Notes:
%       cos(0) = 1 marks the slow oscillation (SO) upstate; cos(+/-pi) = -1
%       marks the SO downstate. Convenience overlay for SO-phase figures.
%
%   See also: plot_stage_prefphase, phasehistogram
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

    xs = -pi:0.01:pi;
    ys = cos(xs);
    hold on;
    plot(xs,ys,'k');
    xlim([-pi pi])
    yline(0,'k--')
    set(gca,'YTick',[],'box','off','YColor','none')
    xticks(-pi:pi/2:pi);
    xticklabels({'-\pi','-\pi/2','0', '\pi/2','\pi'});
    xlabel('SO Phase (rad)','fontsize',10);

    % Add "SO Rising" label
    text(-pi/2-0.4, cos(-pi/2-0.4) + 0.5, 'SO Rising', ...
        'HorizontalAlignment', 'center', 'Rotation', 45, 'FontSize', 10);
    
    % Add "SO Falling" label
    text(pi/2+0.4, cos(pi/2+0.4) + 0.5, 'SO Falling', ...
        'HorizontalAlignment', 'center', 'Rotation', 315, 'FontSize', 10);
    ax1 = gca; 
    ax1.FontSize = 10;

end


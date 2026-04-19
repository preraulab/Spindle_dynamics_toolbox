function [Sp_adj] = FinerModCardinalSpline(ord,c_pt_times_all,s,spline_resol)
%FINERMODCARDINALSPLINE  Modified cardinal spline basis on a finer evaluation grid (non-uniform knots)
%
%   Usage:
%       Sp_adj = FinerModCardinalSpline(ord, c_pt_times_all, s, spline_resol)
%
%   Inputs:
%       ord            : integer - history order (number of base bins) -- required
%       c_pt_times_all : 1xK double - knot locations (bin units) -- required
%       s              : double - cardinal spline tension parameter -- required
%       spline_resol   : double - sub-bin spline evaluation resolution -- required
%
%   Outputs:
%       Sp_adj : (ord/spline_resol)xK double - spline basis matrix evaluated at the finer grid
%
%   Notes:
%       Boundary knots use the reflected tangent trick (Sarmashghi et al.,
%       PLoS One 2021); interior knots use tension-scaled non-uniform
%       spacing factors. Companion to ModifiedCardinalSpline but on a
%       spline_resol sub-bin grid.
%
%   See also: ModifiedCardinalSpline, plot_hist_curve
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%%

lastknot = ord;
numknots = length(c_pt_times_all);
Sp_adj = zeros(lastknot/spline_resol,numknots);

for i= 1:lastknot/spline_resol
   nearest_c_pt_index = max(find(c_pt_times_all<i*spline_resol));
   nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
   next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
   
   u = (i*spline_resol-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
   lb = (c_pt_times_all(3) - c_pt_times_all(1))/(c_pt_times_all(2)-c_pt_times_all(1));
   le = (c_pt_times_all(end) - c_pt_times_all(end-2))/(c_pt_times_all(end) - c_pt_times_all(end-1));   
   
   % Beginning knot 
   if nearest_c_pt_time == c_pt_times_all(1)  % Fixed Mehrad version
           p = [u^3 u^2 u 1]*[2-(s/lb) -2 s/lb;(s/lb)-3 3 -s/lb;0 0 0;1 0 0];
           Sp_adj(i,nearest_c_pt_index:nearest_c_pt_index+2) = p; 
   % End knot
   elseif nearest_c_pt_time==c_pt_times_all(end-1) % Fixed Mehrad version
           p=[u^3 u^2 u 1]*[-s/le 2 -2+(s/le);2*s/le -3 3-(2*s/le);-s/le 0 s/le;0 1 0];
           Sp_adj(i,nearest_c_pt_index-1:nearest_c_pt_index+1) = p;  
   % Interior knots
   else
           prev_c_pt_time = c_pt_times_all(nearest_c_pt_index-1);
           next2 = c_pt_times_all(nearest_c_pt_index+2);
           l1 = (next_c_pt_time-prev_c_pt_time)/(next_c_pt_time-nearest_c_pt_time); % scale factors for non-uniform spacing 
           l2 = (next2-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time); % scale factors for non-uniform spacing
           p=[u^3 u^2 u 1]*[-s/l1 2-s/l2 s/l1-2 s/l2;2*s/l1 s/l2-3 3-2*s/l1 -s/l2;-s/l1 0 s/l1 0;0 1 0 0];
           Sp_adj(i,nearest_c_pt_index-1:nearest_c_pt_index+2) = p;
    end
end


end
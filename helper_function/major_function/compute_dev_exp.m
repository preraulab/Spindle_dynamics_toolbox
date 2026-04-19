function [dev_exp_sta,dev_exp_sop] = compute_dev_exp(BinData)
%COMPUTE_DEV_EXP  Compute proportional deviance explained by each model factor
%
%   Usage:
%       [dev_exp_sta, dev_exp_sop] = compute_dev_exp(BinData)
%
%   Inputs:
%       BinData : struct - all binned data required for the nested model fits -- required
%
%   Outputs:
%       dev_exp_sta : 1x3 double - deviance explained from [stage, phase, history]
%       dev_exp_sop : 1x3 double - deviance explained from [SOP,   phase, history]
%
%   Notes:
%       Fits a sequence of nested Poisson GLMs (null, single-factor, and
%       three-factor) via glmfit and computes the relative contribution of
%       each factor to the total deviance reduction. Accompanies Chen et
%       al., PNAS 2025 (https://doi.org/10.1073/pnas.2405276121).
%
%   See also: specify_mdl, build_design_mt, glmfit
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Fit all component models
% SOphase, stage, SOpower,history
y = BinData.y;

% Baseline, NREM sleep
[ModelSpec0] = specify_mdl([0 0 0 0],[],'hist_choice',[],'control_pt',[]);
X0 = build_design_mt(ModelSpec0,BinData);
[b0, dev0, stats0] = glmfit(X0,y,'poisson'); 
[CIF0,~,~] = glmval(b0,X0,'log',stats0);

% Phase only Model
[ModelSpec1_p] = specify_mdl([1 0 0 0],[],'hist_choice',[],'control_pt',[]);
X1_p = build_design_mt(ModelSpec1_p,BinData);
[b1_p, dev1_p, stats1_p] = glmfit(X1_p,y,'poisson'); 
[CIF1_p,~,~] = glmval(b1_p,X1_p,'log',stats1_p);

% Stage only model
[ModelSpec1_sta] = specify_mdl([0 1 0 0],[],'hist_choice',[],'control_pt',[]);
X1_sta = build_design_mt(ModelSpec1_sta,BinData);
[b1_sta, dev1_sta, stats1_sta] = glmfit(X1_sta,y,'poisson'); 
[CIF1_sta,~,~] = glmval(b1_sta,X1_sta,'log',stats1_sta);

% SOP only model
[ModelSpec1_sop] = specify_mdl([0 0 1 0],[],'hist_choice',[],'control_pt',[]);
X1_sop = build_design_mt(ModelSpec1_sop,BinData);
[b1_sop, dev1_sop, stats1_sop] = glmfit(X1_sop,y,'poisson'); 
[CIF1_sop,~,~] = glmval(b1_sop,X1_sop,'log',stats1_sop);

% History only model
[ModelSpec1_h] = specify_mdl([0 0 0 1],[]);
X1_h = build_design_mt(ModelSpec1_h,BinData);
[b1_h, dev1_h, stats1_h] = glmfit(X1_h,y,'poisson'); 
[CIF1_h,~,~] = glmval(b1_h,X1_h,'log',stats1_h);

% Model with sop, phase and history 
[ModelSpec3_sop] = specify_mdl([1 0 1 1],[]);
X3_sop = build_design_mt(ModelSpec3_sop,BinData);
[b3_sop, dev3_sop, stats3_sop] = glmfit(X3_sop,y,'poisson'); 
[CIF3_sop,~,~] = glmval(b3_sop,X3_sop,'log',stats3_sop);

% Model with stage, phase and history 
[ModelSpec3_sta] = specify_mdl([1 1 0 1],[]);
X3_sta = build_design_mt(ModelSpec3_sta,BinData);
[b3_sta, dev3_sta, stats3_sta] = glmfit(X3_sta,y,'poisson'); 
[CIF3_sta,~,~] = glmval(b3_sta,X3_sta,'log',stats3_sta);

%% Compute relative contribution (proportional dev explained) for each factor
% Deviance explained can be derived from models with stage
devr_all_sta = dev0 - dev3_sta;     % Total dev reduction from null model
devr_sta = dev0 - dev1_sta;         % Dev reduction from null to model with stage
devr_h = dev0 - dev1_h;
devr_p = dev0 - dev1_p;
dev_exp_sta = [devr_sta devr_p devr_h]./devr_all_sta;

% Or from models with SOP
devr_all_sop = dev0 - dev3_sop;
devr_sop = dev0 - dev1_sop;
devr_h = dev0 - dev1_h;
devr_p = dev0 - dev1_p;
dev_exp_sop = [devr_sop devr_p devr_h]./devr_all_sop;

end
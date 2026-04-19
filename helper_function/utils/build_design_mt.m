function [X] = build_design_mt(ModelSpec,BinData)
%BUILD_DESIGN_MT  Construct a GLM design matrix from user-selected factors and interactions
%
%   Usage:
%       X = build_design_mt(ModelSpec, BinData)
%
%   Inputs:
%       ModelSpec : struct - model specification (BinarySelect, InteractSelect, ...) -- required
%       BinData   : struct - all binned data (sta, sop, phase, sp_hist, RW) -- required
%
%   Outputs:
%       X : TxP double - design matrix (size depends on selected factors and interactions)
%
%   Notes:
%       Factors have fixed order [SOphase, stage, SOpower, history]. SOP is
%       rescaled internally by /5 and quadratic terms by /10 to improve
%       numerical conditioning. Supported interactions: stage:SOphase,
%       SOpower:SOphase, and combined stage:SOphase + stage:history.
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: specify_mdl, preprocessToDesignMatrix, ismember_interaction
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% Deal with inputs

% Extract model specifications
interactions = ModelSpec.InteractSelect;
BinarySelectStr = num2str(ModelSpec.BinarySelect, '%d');
BinarySelectStr = strrep(BinarySelectStr, ' ', ''); % Convert binary array to a string number

% Extract factors
sta = BinData.sta;
sop = BinData.sop;
phase = BinData.phase;
sp_hist = BinData.sp_hist;
RW = BinData.RW;

% Rescale SOP to avoid numerical issues in model fitting
% Alternative: Standardize SOP, more robust 
sop = sop/5;
scale_factor = 10;

% Initialize X
X = [];
%% Construct design matrix
% 1) Without interactions______________________________________________________
if isempty(interactions)
    switch BinarySelectStr
        % Single factor
        case '0000' % No factors, NREM sleep
            X = RW;                             
        case '1000' % phase alone
            X = [RW cos(phase) sin(phase)];     
        case '0100' % stage alone
            X = sta(:,[1 3:5]);                
        case '0010' % sop alone
            X = [RW sop (sop.^2)/scale_factor]; 
        case '0001' % history alone
            X = [RW sp_hist];                   
        % Two-factor
        case '1100' % phase + stage
            X = [sta(:,[1 3:5]) cos(phase) sin(phase)]; 
        case '1010' % phase + sop
            X = [RW sop (sop.^2)/scale_factor cos(phase) sin(phase)]; 
        case '1001' % phase + hist
            X = [RW cos(phase) sin(phase) sp_hist];      
        case '0101' % stage + hist
            X = [sta(:,[1 3:5]) sp_hist];                  
        case '0011' % sop + hist
            X = [RW sop (sop.^2)/scale_factor sp_hist];  
        % Three-factor
        case '1101' % sta,phase,hist
            X = [sta(:,[1 3:5]) cos(phase) sin(phase) sp_hist];
        case '1011' % sop,phase,hist
            X = [RW sop (sop.^2)/scale_factor cos(phase) sin(phase) sp_hist];
        otherwise
            error('The model components are not properly specified.')
    end

% 2) With only stage:SOphase interaction ______________________________________________________
elseif all(ismember_interaction(interactions,{'stage:SOphase'}))
     switch BinarySelectStr
        case '1100' % phase + stage
            X = [cos(phase) sin(phase) sta(:,[1 3:5]) sta(:,[1 3:5]).*cos(phase) sta(:,[1 3:5]).*sin(phase)];
        case '1101' % sta,phase,hist
            X = [cos(phase) sin(phase) sta(:,[1 3:5]) sta(:,[1 3:5]).*cos(phase) sta(:,[1 3:5]).*sin(phase) sp_hist];
     end

% 3) With only SOpower:SOphase interaction ______________________________________________________
elseif all(ismember_interaction(interactions,{'SOpower:SOphase'}))
    switch BinarySelectStr
        case '1010' % phase + sop
            X = [RW sop (sop.^2)/scale_factor cos(phase) sin(phase) ...
                 cos(phase).*sop  cos(phase).*(sop.^2)/scale_factor ...
                 sin(phase).*sop  sin(phase).*(sop.^2)/scale_factor];
        case '1011' % sop,phase,hist
            X = [RW sop (sop.^2)/scale_factor cos(phase) sin(phase) ...
                 cos(phase).*sop  cos(phase).*(sop.^2)/scale_factor ...
                 sin(phase).*sop  sin(phase).*(sop.^2)/scale_factor sp_hist];
    end

 % 4) With both stage：SOphase AND stage：history interactions ___________________________
 elseif all(ismember_interaction(interactions,{'stage:SOphase','stage:history'}))&&all(ModelSpec.BinarySelect == [1 1 0 1])
            X = [cos(phase) sin(phase) sta(:,[1 3:5]) sta(:,[1 3:5]).*cos(phase) sta(:,[1 3:5]).*sin(phase) ...
                 sp_hist sta(:,2).*sp_hist sta(:,3).*sp_hist];
 else
        error('Model component or interaction term is not properly specified')
end
 
if isempty(X)
    error('Model component or interaction term is not properly specified')
end

end

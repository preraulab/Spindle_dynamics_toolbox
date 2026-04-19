%PHASEHISTOGRAM  Draw a polar histogram with mean-vector arrow for phase data
%
%   Usage:
%       [theta_mean, rho_mean, h_phist, h_pax, h_ml] = phasehistogram(phases, amps, 'Name', Value, ...)
%
%   Inputs:
%       phases : 1xN double - phase values (rad) -- required
%       amps   : 1xN double - amplitudes used to weight the mean vector (default: ones)
%
%   Name-Value Pairs:
%       <polarhistogram options> : any name-value accepted by polarhistogram() (NumBins, FaceColor, FaceAlpha, ...)
%
%   Outputs:
%       theta_mean : double - mean angle (rad)
%       rho_mean   : double - mean magnitude (in histogram units after normalization to 'pdf')
%       h_phist    : handle - polarhistogram handle
%       h_pax      : handle - polar axes handle
%       h_ml       : handle - mean-line handle
%
%   Example:
%       phases = mod(randn(1,1000) + pi/2, 2*pi);
%       figure; phasehistogram(phases, 1, 'NumBins', 25, 'FaceColor', 'blue', 'FaceAlpha', 0.3);
%
%   See also: polarhistogram, polaraxes
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

function [theta_mean, rho_mean, h_phist, h_pax, h_ml] = phasehistogram(phases, amps, varargin)
if nargin==0
    error('Must input phases');
end

if nargin<2 || isempty(amps)
    amps = ones(size(phases));
end

%Compute the mean population vector
vect_mean = mean(amps.*exp(1i*phases));

%Get the mean magnitude and angle
rho_mean = abs(vect_mean);
theta_mean = angle(vect_mean);

%Set default to normalization
varargin = [{'Normalization'}, {'pdf'}, varargin(:)'];

%Plot histogram
h_phist = polarhistogram(phases,varargin{:});
h_pax = gca;
h_pax.ThetaAxisUnits = 'radians';
%h_pax.ThetaTick = 0:pi/4:2*pi;
%h_pax.ThetaTick = 0:pi/2:2*pi;
h_pax.ThetaTick = 0:pi:2*pi;
%h_pax.ThetaTickLabel = {'0','\pi/4','\pi/2','3\pi/4' '\pm\pi','-3\pi/4', '-\pi/2','-\pi/4'};
%h_pax.ThetaTickLabel = {'0','\pi/2', '\pm\pi', '-\pi/2'};
h_pax.ThetaTickLabel = {'0','\pm\pi'};
h_pax.FontSize = 12;

%Add mean arrow
hold on
h_ml = polarplot([theta_mean theta_mean],[0 rho_mean],'linestyle','-','color','r','linewidth',2);


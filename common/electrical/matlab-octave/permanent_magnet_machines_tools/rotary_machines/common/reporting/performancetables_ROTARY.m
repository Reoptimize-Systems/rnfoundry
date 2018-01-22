function [ptables, pdata] = performancetables_ROTARY(design, simoptions, rpm, LoadVals, outfields, varargin)
% generates tables of performance data a multiple speed and load points for
% a rotary machine design
%
% Syntax
%
% [ptables] = performancetables_ROTARY(design, simoptions, rpm, LoadVals, outfields)
%
% Input
%
%   design, simoptions - the machine design and simoptions structures
%
%   rpm - vector of rpm values
%
%   LoadVals - vector of load resistance to phase resistance ratios
%
%   outfields - optional cell array of strings containing field names from
%     the simulation output design array which will be added to the same
%     named field in the output table structure.
%
% Output
%
%   ptables - structure containing tables of output data at each
%     corresponding combination of speed and loading. The fields of ptables
%     will contain at least the following:
%
%       OmegaPeak
%       TorquePtoMean
%       TorquePtoPeak
%       PowerLoadMean
%       Efficiency
%       IPhasePeak
%       IPhaseRms
%       ICoilPeak
%       ICoilRms
%       EMFPhasePeak
%       EMFPhaseRms
%       JCoilPeak
%       JCoilRms
%       EnergyLoadTotal
%       PowerPhaseRMean
%       TorquePtoMean
%       TorquePtoPeak
%       PowerLossEddyMean
%       PowerLossMean
%       FrequencyPeak
%
%     each table will be a (n x m) matrix where n is the number of speed
%     points and m is the number of load ratio points supplied.
%

% Created by Richard Crozier 2013-2015

    options.UseParFor = false;
    options.LoadSpecType = 'ratio';
    options.DoSimFun = true;
    options.Verbose = false;
    
    options = parse_pv_pairs (options, varargin);
    
    if nargin < 5   
        outfields = {};
    end

	% append to any machine specific output fields already passed in
    outfields = [ outfields, { 'OmegaPeak', ...
                               'TorquePtoMean', ...
                               'TorquePtoPeak' } ];

    [ptables, pdata] = performancetables_AM ( design, simoptions, rpm, LoadVals, outfields, ...
                    'UseParFor', options.UseParFor, ...
                    'LoadSpecType', options.LoadSpecType, ...
                    'SpeedSetupFcn', @setrpm, ...
                    'DoSimFun', options.DoSimFun, ...
                    'Verbose', options.Verbose);

end

function [design, simoptions] = setrpm (design, simoptions, rpm, options)
    simoptions.RPM = rpm;
end

function ptables = performancetables_RADIAL(design, simoptions, rpm, RgVRc, outfields)
% generates tables of performance data a multiple speed and load points for
% a radial flux rotary machine design
%
% Syntax
%
% ptables = performancetables_RADIAL(design, simoptions, rpm, RgVRc, outfields)
%
% Input
%
%   design, simoptions - the machine design and simoptions structures
%
%   rpm - vector of rpm values
%
%   RgVRc - vector of load resistance to phase resistance ratios
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
%       MaxTorquePto
%       PowerLossEddyMean
%       PowerLossMean
%       FrequencyPeak
%
%     each table will be a (n x m) matrix where n is the number of speed
%     points and m is the number of load ratio points supplied.
%

% Created by Richard Crozier 2013

    if nargin < 5
        outfields = {};
    end

    [ptables] = performancetables_ROTARY(design, simoptions, rpm, RgVRc, outfields);

end
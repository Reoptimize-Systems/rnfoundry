function [ptables] = performancetables_ROTARY(design, simoptions, rpm, RlVRp, outfields)
% generates tables of performance data a multiple speed and load points for
% a rotary machine design
%
% Syntax
%
% [ptables] = performancetables_ROTARY(design, simoptions, rpm, RlVRp, outfields)
%
% Input
%
%   design, simoptions - the machine design and simoptions structures
%
%   rpm - vector of rpm values
%
%   RlVRp - vector of load resistance to phase resistance ratios
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

    radialfields = { 'PowerLoadMean', ...
                     'Efficiency', ...
                     'IPhasePeak', ...
                     'IPhaseRms', ...
                     'ICoilPeak', ...
                     'ICoilRms', ...
                     'EMFPhasePeak', ...
                     'EMFPhaseRms', ...
                     'JCoilPeak', ...
                     'JCoilRms', ...
                     'EnergyLoadTotal', ...
                     'PowerPhaseRMean', ...
                     'MaxTorquePto', ...
                     'PowerLossEddyMean', ...
                     'PowerLossMean', ...
                     'FrequencyPeak' };
                
    outfields = [ outfields, radialfields ];
    
    % preallocate the fields
    
    isfirstrun = true;
    
    for rpmind = 1:numel(rpm)
        
        for RgVRcind  = 1:numel(RlVRp)
            
            simoptions.RPM = rpm(rpmind);
            simoptions.RlVRp = RlVRp(RgVRcind);
            design.RlVRp = RlVRp(RgVRcind);
            simoptions.abstol = [];
            
            % simulate the machine
            [~, ~, ~, design, simoptions] = simulatemachine_AM( design, ...
                                                                simoptions, ...
                                                                simoptions.simfun, ...
                                                                simoptions.finfun, ...
                                                                simoptions.odeevfun, ...
                                                                simoptions.resfun );
                                                                 
            if isfirstrun
                isfirstrun = false;
                % remove simfun, so we don't repeat fea etc. on subsequent
                % runs
                simoptions.simfun = [];
            end
            
            % copy the output over to the performance tables
            for find = 1:numel(outfields)
                ptables.(outfields{find})(rpmind,RgVRcind) = design.(outfields{find});
            end
            
        end
        
    end

end
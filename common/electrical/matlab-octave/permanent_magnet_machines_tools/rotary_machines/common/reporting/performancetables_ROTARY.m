function [ptables, pdata] = performancetables_ROTARY(design, simoptions, rpm, RlVRp, outfields)
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

    if nargin < 5
        outfields = {};
    end
    
    % get output fields common to all rotary machines
    rotaryfields = { 'PowerLoadMean', ...
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
                     'TorquePtoMean', ...
                     'TorquePtoPeak', ...
                     'PowerLossEddyMean', ...
                     'PowerLossMean', ...
                     'FrequencyPeak' };
                
	% append to any machine specific output fields already passed in
    outfields = [ outfields, rotaryfields ];

    % do design data gathering function if necessary
    if ~isempty(simoptions.simfun)
        % Analyse the machine and gather desired data
        [design, simoptions] = feval(simoptions.simfun, design, simoptions);
    end
    
	% remove simfun, so we don't repeat fea etc. on subsequent runs
	simoptions.simfun = [];

    % preallocate the fields
    ptables = struct ();
    pdata = nan * ones (numel(rpm) * numel (RlVRp), numel (outfields));
    rpmslice = nan * ones (numel(rpm) * numel (RlVRp), 1);
    RlVRpslice = rpmslice;
    
    rpmRlVRpind = 1;
    for rpmind = 1:numel(rpm)
        
        for RgVRcind = 1:numel(RlVRp)
            
            rpmslice(rpmRlVRpind) = rpm(rpmind);
            RlVRpslice(rpmRlVRpind) = RlVRp(RgVRcind);
            
            rpmRlVRpind = rpmRlVRpind + 1;
            
        end
        
    end
    
    % perform the simulations and gather the data
    parfor rpmRlVRpind = 1:rpmRlVRpind-1
        pdata(rpmRlVRpind,:) = simfcn (design, simoptions, outfields, rpmslice(rpmRlVRpind), RlVRpslice(rpmRlVRpind));
    end
    
    % construct the output tables from the matrix of data
    for rpmind = 1:numel(rpm)
        
        for RgVRcind = 1:numel(RlVRp)
            
            % copy the output over to the performance tables
            for find = 1:numel(outfields)
                ptables.(outfields{find})(rpmind,RgVRcind) = pdata((rpmind-1)*numel(RlVRp) + RgVRcind, find);
            end
            
        end
        
    end

end

function pdata = simfcn (design, simoptions, outfields, rpm, RlVRp)

        simoptions.RPM = rpm;
        simoptions.RlVRp = RlVRp;
        design.RlVRp = RlVRp;
        design = rmiffield (design, 'LoadResistance');
        design = rmiffield (design, 'slm_fluxlinkage');
        simoptions.abstol = [];
        
        % simulate the machine
        [~, ~, ~, design, ~] = simulatemachine_AM( design, ...
                                                   simoptions, ...
                                                   simoptions.simfun, ...
                                                   simoptions.finfun, ...
                                                   simoptions.odeevfun, ...
                                                   simoptions.resfun );

        % copy the output over to the performance tables
        pdata = nan * ones (1,numel(outfields));
        for find = 1:numel(outfields)
            pdata(find) = design.(outfields{find})(1);
        end
            
end
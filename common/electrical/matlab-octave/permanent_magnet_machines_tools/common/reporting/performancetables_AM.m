function [ptables, pdata] = performancetables_AM (design, simoptions, speeds, LoadVals, outfields, varargin)
% generates tables of performance data a multiple speed and load points for
% a machine design
%
% Syntax
%
% [ptables] = performancetables_AM (design, simoptions, speeds, LoadVals, outfields)
% [ptables] = performancetables_AM (..., 'Parameter', Value)
%
% Input
%
%   design, simoptions - the machine design and simoptions structures
%
%   speeds - vector of speed values
%
%   LoadVals - vector of load resistances. Be defualt these are specified
%     as load to phase resistance ratios. 
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

    options.UseParFor = false;
    options.LoadSpecType = 'ratio';
    options.SpeedSetupFcn = @setrpm;
    options.DoSimFun = true;
    
    options = parse_pv_pairs (options, varargin);
    
    if nargin < 5
        outfields = {};
    end
    
    % get output fields common to all rotary machines
    commonfields = { 'PowerLoadMean', ...
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
                     'FrequencyPeak', ...
                     'PowerFactorEstimate' };
                
	% append to any machine specific output fields already passed in
    outfields = [ outfields, commonfields ];

    % do design data gathering function if necessary
    if options.DoSimFun && ~isempty(simoptions.simfun)
        % Analyse the machine and gather desired data
        [design, simoptions] = feval (simoptions.simfun, design, simoptions);
    end
    
	% remove simfun, so we don't repeat fea etc. on subsequent runs
	simoptions.simfun = [];

    % preallocate the fields
    ptables = struct ();
    pdata = nan * ones (numel(speeds) * numel (LoadVals), numel (outfields));
    rpmslice = nan * ones (numel(speeds) * numel (LoadVals), 1);
    LoadValsslice = rpmslice;
    
    rpmLoadValsind = 1;
    for rpmind = 1:numel(speeds)
        
        for LoadValsind = 1:numel(LoadVals)
            
            rpmslice(rpmLoadValsind) = speeds(rpmind);
            LoadValsslice(rpmLoadValsind) = LoadVals(LoadValsind);
            
            rpmLoadValsind = rpmLoadValsind + 1;
            
        end
        
    end
    
    % perform the simulations and gather the data
    if options.UseParFor
        parfor rpmLoadValsind = 1:rpmLoadValsind-1
            pdata(rpmLoadValsind,:) = simfcn (design, simoptions, outfields, rpmslice(rpmLoadValsind), LoadValsslice(rpmLoadValsind), options);
        end
    else
        for rpmLoadValsind = 1:rpmLoadValsind-1
            pdata(rpmLoadValsind,:) = simfcn (design, simoptions, outfields, rpmslice(rpmLoadValsind), LoadValsslice(rpmLoadValsind), options);
        end
    end
    
    % construct the output tables from the matrix of data
    for rpmind = 1:numel(speeds)
        
        for LoadValsind = 1:numel(LoadVals)
            
            % copy the output over to the performance tables
            for find = 1:numel(outfields)
                if isfield (design, outfields{find})
                    ptables.(outfields{find})(rpmind,LoadValsind) = pdata((rpmind-1)*numel(LoadVals) + LoadValsind, find);
                end
            end
            
        end
        
    end
    
    % report a bit of information about what was done
    ptables.SimulationInfo.Speeds = speeds;
    ptables.SimulationInfo.LoadVals = LoadVals;
    ptables.SimulationInfo.LoadSpecType = options.LoadSpecType;
    if isa (options.SpeedSetupFcn, 'function_handle')
        ptables.SimulationInfo.SpeedSetupFcn = func2str (options.SpeedSetupFcn);
    else
        ptables.SimulationInfo.SpeedSetupFcn = options.SpeedSetupFcn;
    end

end

function [design, simoptions] = setrpm (design, simoptions, rpm, options)
    simoptions.RPM = rpm;
end

function pdata = simfcn (design, simoptions, outfields, rpm, LoadVal, options)

    [design, simoptions] = feval (options.SpeedSetupFcn, design, simoptions, rpm, options);

    design = rmiffield (design, 'LoadResistance');
    design = rmiffield (design, 'RlVRp');

    switch options.LoadSpecType
        case 'ratio'
            simoptions.RlVRp = LoadVal;
            design.RlVRp = LoadVal;
        case 'resistance'
            simoptions.LoadResistance = LoadVal;
            design.LoadResistance = LoadVal;
        otherwise
            error ('Unrecognised specification type')
    end

    % remove flux linkage slm to trigger rerun of postprocessing
    % function
    design = rmiffield (design, 'slm_fluxlinkage');
    simoptions.abstol = [];

    % simulate the machine
    [~, ~, ~, design, ~] = simulatemachine_AM ( design, ...
                                                simoptions );

    % copy the output over to the performance tables
    pdata = nan * ones (1,numel(outfields));
    for find = 1:numel(outfields)
        if isfield (design, outfields{find})
            pdata(find) = design.(outfields{find})(1);
        end
    end  
end
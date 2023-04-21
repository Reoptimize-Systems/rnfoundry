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

    opts.UseParFor = false;
    opts.LoadSpecType = 'ratio';
    opts.SpeedSetupFcn = @setrpm;
    opts.DoSimFun = true;
    opts.Verbose = false;
    
    opts = parse_pv_pairs (opts, varargin);
    
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
    if opts.DoSimFun && ~isempty(simoptions.ODESim.PreProcFcn)
        % Analyse the machine and gather desired data
        [design, simoptions] = feval (simoptions.ODESim.PreProcFcn, design, simoptions);
    end
    
	% remove simfun, so we don't repeat fea etc. on subsequent runs
	simoptions.ODESim.PreProcFcn = [];

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
    nsims = rpmLoadValsind - 1;
    
    % report a bit of information about what was done
    ptables.SimulationInfo.Speeds = speeds;
    ptables.SimulationInfo.LoadVals = LoadVals;
    ptables.SimulationInfo.LoadSpecType = opts.LoadSpecType;
    if isa (opts.SpeedSetupFcn, 'function_handle')
        ptables.SimulationInfo.SpeedSetupFcn = func2str (opts.SpeedSetupFcn);
    else
        ptables.SimulationInfo.SpeedSetupFcn = opts.SpeedSetupFcn;
    end
    ptables.SimulationInfo.Design = design;
    ptables.SimulationInfo.Simoptions = simoptions;
    
    % perform the simulations and gather the data
    if opts.UseParFor
        parfor rpmLoadValsind = 1:nsims-1
            if opts.Verbose, fprintf (1, 'Performing sim %d of %d', rpmLoadValsind, nsims); end
            pdata(rpmLoadValsind,:) = simfcn (design, simoptions, outfields, rpmslice(rpmLoadValsind), LoadValsslice(rpmLoadValsind), opts);
        end
    else
        % do all but last sim, which is done later (getting desing
        % stucture)
        for rpmLoadValsind = 1:nsims-1
            if opts.Verbose, fprintf (1, 'Performing sim %d of %d', rpmLoadValsind, nsims); end
            pdata(rpmLoadValsind,:) = simfcn (design, simoptions, outfields, rpmslice(rpmLoadValsind), LoadValsslice(rpmLoadValsind), opts);
        end
    end
    
    % get last data, and get design containing output fields
    if opts.Verbose, fprintf (1, 'Performing sim %d of %d', nsims, nsims); end
    [pdata(nsims,:), design] = simfcn (design, simoptions, outfields, rpmslice(nsims), LoadValsslice(nsims), opts);
    
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

end

function [design, simoptions] = setrpm (design, simoptions, rpm, options)
    simoptions.RPM = rpm;
end

function [pdata, design] = simfcn (design, simoptions, outfields, rpm, LoadVal, options)

    [design, simoptions] = feval (options.SpeedSetupFcn, design, simoptions, rpm, options);

    design = rmiffield (design, 'LoadResistance');
    design = rmiffield (design, 'RlVRp');

    % set PostPreProcessingComplete = false so that the load is
    % recalculated properly if necessary
    design.PostPreProcessingComplete = false;

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
%     design = rmiffield (design, 'slm_fluxlinkage');
    simoptions.ODESim.AbsTol = [];

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
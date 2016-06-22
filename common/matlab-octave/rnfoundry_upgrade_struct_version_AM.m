function [design, simoptions] = rnfoundry_upgrade_struct_version_AM (design, simoptions)
% this version converts design and simoptions structures to use the new
% field names used for various options
%

    subfname = 'ODESim';

    fnames = { 'simfun', 'finfun', 'resfun', 'odeevfun', ...
               'torquefcn', 'torquefcnargs', 'forcefcn', ...
               'forcefcnargs', 'maxstep', 'abstol', ...
               'reltol', 'skip' };
           
    newfnames = { 'PreProcFcn', 'PostPreProcFcn', 'PostSimFcn', 'EvalFcn', ...
                  'TorqueFcn', 'TorqueFcnArgs', 'ForceFcn', ...
                  'ForceFcnArgs', 'MaxStep', 'AbsTol', ...
                  'RelTol', 'Skip' };
    
    simoptions = upgradefield (simoptions, fnames, newfnames, subfname);
    
    
    subfname = 'BuoySim';

    fnames = { 'Buoy', 'BuoyParameters', 'SeaParameters', ...
               'buoylibdir', 'HeaveFile', 'SurgeFile', ...
               'ExcitationFile', 'HydroCoeffsFile', ...
               'NRadiationCoefs', 'EndStopks', 'maxAllowedxT', ...
               'tether_length'};
           
    newfnames = fnames;
    
    simoptions = upgradefield (simoptions, fnames, newfnames, subfname);

end

function S = upgradefield (S, fnames, newfnames, subfname)

    for ind = 1:numel (fnames)
        
        if nargin > 3

            S = setfieldifabsent (S, subfname, struct ());

            if isfield (S, fnames{ind})
                S.(subfname).(newfnames{ind}) = S.(fnames{ind});
                S = rmfield (S, fnames{ind});
            end

        else

            if isfield (S, fnames{ind})
                S.(newfnames{ind}) = S.(fnames{ind});
                S = rmfield (S, fnames{ind});
            end

        end
    
    end
    
end
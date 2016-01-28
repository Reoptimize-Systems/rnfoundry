function namestr = optname_AM (simoptions, startstr)
% creates a descriptive name string for common optimisation aspects
    
    if nargin < 2
        startstr = '';
    end
    
    namestr = startstr;
    
    if isfield (simoptions, 'TargetPowerLoadMean')
        
        [num, prefix] = scalar_field_num_and_si_prefix (simoptions, 'TargetPowerLoadMean');
        
        namestr = [ namestr, '_' num, prefix, 'W' ];
    end
    
    if isfield (simoptions, 'min_EMFPhaseRms')
        
        [num, prefix] = scalar_field_num_and_si_prefix (simoptions, 'min_EMFPhaseRms');
        
        namestr = [ namestr, '_+' num, prefix, 'V' ];
    end
    
    if isfield (simoptions, 'max_EMFPhaseRms')
        
        [num, prefix] = scalar_field_num_and_si_prefix (simoptions, 'max_EMFPhaseRms');
        
        namestr = [ namestr, '_-' num, prefix, 'V' ];
    end
    
    if isfield (simoptions, 'max_JCoilRms')
        
        [num, prefix] = scalar_field_num_and_si_prefix (simoptions, 'max_JCoilRms');
        
        namestr = [ namestr, '_-' num, prefix, 'Apm2' ];
    end
    
    if isfield (simoptions, 'RlVRp')
        
        namestr = [ namestr, '_RlVRp_' num2str(simoptions.RlVRp) ];
    end
    
    namestr = strrep (namestr, '.', 'pt');

end
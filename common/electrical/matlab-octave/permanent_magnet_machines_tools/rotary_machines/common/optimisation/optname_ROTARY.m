function namestr = optname_ROTARY (simoptions, startstr)
% creates a descriptive name string for optimisation aspects common to
% rotary machines

    if nargin < 2
        startstr = '';
    end
    
    namestr = startstr;
    
    if isfield (simoptions, 'RPM')
        namestr = [ namestr, '_' num2str(simoptions.RPM), 'rpm' ];
    end
    
    namestr = optname_AM (simoptions, namestr);

end
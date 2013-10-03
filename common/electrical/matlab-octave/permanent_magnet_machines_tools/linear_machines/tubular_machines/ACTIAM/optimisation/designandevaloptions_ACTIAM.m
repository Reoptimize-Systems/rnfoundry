function evaloptions = designandevaloptions_ACTIAM(evaloptions)

    if nargin == 0
        evaloptions = [];
    end
    
    evaloptions = designandevaloptions_linear(evaloptions);
    
    if isempty(evaloptions)
        
            evaloptions.E = [207e9 100e6 207e9];
            evaloptions.StructSections = 10;
            evaloptions.vCoil = 0.3;
            evaloptions.vSheath = 0.28;
            evaloptions.coilYieldStrength = 70e6;
            evaloptions.presimfinfun = 'prescribedmotfinfun_ACTIAM';

    else
            % If the evaloptions structure is supplied, we check each possible
            % optional values and set any not supplied to the default value

            if ~isfield(evaloptions, 'E')
                evaloptions.E = [207e9 100e6 207e9];
            end

            if ~isfield(evaloptions, 'StructSections')
                evaloptions.StructSections = 10;
            end

            if ~isfield(evaloptions, 'vCoil')
                evaloptions.vCoil = 0.3;
            end

            if ~isfield(evaloptions, 'vSheath')
                evaloptions.vSheath = 0.28;
            end

            if ~isfield(evaloptions, 'coilYieldStrength')
                evaloptions.coilYieldStrength = 70e6;
            end  

            evaloptions.presimfinfun = 'prescribedmotfinfun_ACTIAM';
            
    end
end
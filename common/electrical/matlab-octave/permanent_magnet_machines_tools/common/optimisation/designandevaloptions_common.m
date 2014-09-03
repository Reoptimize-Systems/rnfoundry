function evaloptions = designandevaloptions_common(evaloptions)
% sets up common default options for electrical machine evaluations

% Copyright Richard Crozer, The University of Edinburgh

    if nargin == 0
        evaloptions = [];
    end
    
    % make an empty struct if evaloptions is empty
    if isempty(evaloptions)
        
        evaloptions = struct;
        
    end
    
    % If the evaloptions structure is supplied, we check each possible
    % optional values and set any not supplied to the default value
    evaloptions = setfieldifabsent(evaloptions, 'MagnetCost', 65);
    
    evaloptions = setfieldifabsent(evaloptions, 'CopperCost', 10);
    
    evaloptions = setfieldifabsent(evaloptions, 'FieldIronCost', 5);
    
    evaloptions = setfieldifabsent(evaloptions, 'ArmatureIronCost', 3);
    
    evaloptions = setfieldifabsent(evaloptions, 'StructMaterialCost', 3);
    
    evaloptions = setfieldifabsent(evaloptions, 'BuoyMassCost', 0.2);

    evaloptions = setfieldifabsent(evaloptions, 'BuoyMassDensity', 1100);

    evaloptions = setfieldifabsent(evaloptions, 'EpoxyCost', 0);

    evaloptions = setfieldifabsent(evaloptions, 'CopperDensity', 8960);
    
    evaloptions = setfieldifabsent(evaloptions, 'ArmatureIronDensity', 7800);
    
    evaloptions = setfieldifabsent(evaloptions, 'FieldIronDensity', 7800);
    
    evaloptions = setfieldifabsent(evaloptions, 'MagnetDensity', 7500);
    
    evaloptions = setfieldifabsent(evaloptions, 'StructMaterialDensity', 7500);
    
    evaloptions = setfieldifabsent(evaloptions, 'StructModulusOfElasticity', 207e9);
    
    evaloptions = setfieldifabsent(evaloptions, 'StructPoissonRatio', 0.31);
    
    evaloptions = setfieldifabsent(evaloptions, 'EpoxyDensity', 2090);
    
    evaloptions = setfieldifabsent(evaloptions, 'CapacityFactor', 0.4);
    
    evaloptions = setfieldifabsent(evaloptions, 'ProjectYears', 30);
    
    evaloptions = setfieldifabsent(evaloptions, 'DiscountRate', 0.08);
    
    evaloptions = setfieldifabsent(evaloptions, 'CostScaleFactor', 100);
    
    evaloptions = setfieldifabsent(evaloptions, 'SkipStructural', false);
    
end

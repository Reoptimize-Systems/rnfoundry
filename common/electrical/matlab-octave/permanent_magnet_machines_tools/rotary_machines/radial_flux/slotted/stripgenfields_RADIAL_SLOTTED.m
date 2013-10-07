function design = stripgenfields_RADIAL_SLOTTED(design)
    
    design = rmiffield(design, 'RsoVRmi');
    design = rmiffield(design, 'RstbVRso');
    design = rmiffield(design, 'tsgVtsb');
    design = rmiffield(design, 'RyoVRstb');
    design = rmiffield(design, 'RyiVRyo');
    design = rmiffield(design, 'RtsgVRso');
    design = rmiffield(design, 'thetamVthetap');
    design = rmiffield(design, 'thetacVthetas');
    design = rmiffield(design, 'thetasgVthetac');
    
    design = stripgenfields_RADIAL(design);
    
end

function design = stripgenfields_RADIAL(design)
    
    design = rmiffield(design, 'lsVtm');
    design = rmiffield(design, 'RmoVRbo');
    design = rmiffield(design, 'RmiVRmo');
    
    design = stripgenfields_ROTARY(design);
    
end

function design = stripgenfields_ROTARY(design)
    
    design = rmiffield(design, 'MaxTorquePto');
    design = rmiffield(design, 'FirstSlotCenter');
    design = rmiffield(design, 'outermagsep');
    
    design = stripgenerated_AM(design);
    
end

function design = stripgenerated_AM(design)
    
    design = rmiffield(design, 'Efficiency');
    design = rmiffield(design, 'BaseScore');
    design = rmiffield(design, 'CoilResistance');
    design = rmiffield(design, 'CoilInductance');
    design = rmiffield(design, 'PhaseResistance');
    design = rmiffield(design, 'PhaseInductance');
    design = rmiffield(design, 'ConductorArea');
    design = rmiffield(design, 'CoilTurns');
    design = rmiffield(design, 'Ntot');
    design = rmiffield(design, 'FEMMTol');
    design = rmiffield(design, 'PoleWidth');
    design = rmiffield(design, 'GridInductance');
    design = rmiffield(design, 'GridResistance');
    design = rmiffield(design, 'AlphaResistivity');
    design = rmiffield(design, 'BranchFac');
    design = rmiffield(design, 'PowerPoles');
    design = rmiffield(design, 'IPhasePeak');
    design = rmiffield(design, 'IPhaseRms');
    design = rmiffield(design, 'EMFPhaseRms');
    design = rmiffield(design, 'JCoilRms');
    design = rmiffield(design, 'JCoilPeak');
    design = rmiffield(design, 'EMFPhasePeak');
    design = rmiffield(design, 'PowerLoadMean');
    design = rmiffield(design, 'PowerPhaseRMean');
    design = rmiffield(design, 'PowerInputMean');
    design = rmiffield(design, 'Efficiency');
    design = rmiffield(design, 'MaxFpto');
    design = rmiffield(design, 'DcAreaFac');
    design = rmiffield(design, 'CostEstimate');
    design = rmiffield(design, 'FemmProblem');
    design = rmiffield(design, 'psilookup');
    design = rmiffield(design, 'CoreLossSLMs');
    design = rmiffield(design, 'slm_eddysfdpart');
    design = rmiffield(design, 'RDCPhase');
    design = rmiffield(design, 'RLoad');
    design = rmiffield(design, 'L');
    design = rmiffield(design, 'CoilPositions');
    design = rmiffield(design, 'slm_psidot');
    design = rmiffield(design, 'slm_fluxlinkage');
    design = rmiffield(design, 'p_gforce');
    design = rmiffield(design, 'gforce');
    design = rmiffield(design, 'gvar');
    design = rmiffield(design, 'coillabellocs');
    design = rmiffield(design, 'yokenodeids');
    design = rmiffield(design, 'feapos');
    design = rmiffield(design, 'ICoilRms');
    
    design = rmiffield(design, 'CoilArea');
    design = rmiffield(design, 'ArmatureIronAreaPerPole');
    design = rmiffield(design, 'MTL');
    design = rmiffield(design, 'Maxdlambdadx');
    design = rmiffield(design, 'EnergyLoadTotal');
    design = rmiffield(design, 'EnergyLoadMean');
    design = rmiffield(design, 'PowerSystemMean');
    design = rmiffield(design, 'PowerLoadPeak');
    design = rmiffield(design, '');
    design = rmiffield(design, '');
    design = rmiffield(design, '');
    
    
    fnames = fieldnames(design);
    
    for i = 1:numel(fnames)
        
        if strncmpi(fliplr(fnames{i}), 'ytlanep', 7)
            design = rmiffield(design, fnames{i});
        end
        
    end
    
    for i = 1:numel(fnames)
        
        if strncmpi(fliplr(fnames{i}), 'tsoC', 4)
            design = rmiffield(design, fnames{i});
        end
        
    end
    
    for i = 1:numel(fnames)
        
        if strncmpi(fliplr(fnames{i}), 'ssaM', 4)
            design = rmiffield(design, fnames{i});
        end
        
    end
    
    for i = 1:numel(fnames)
        
        if strncmpi(fliplr(fnames{i}), 'emuloV', 6)
            design = rmiffield(design, fnames{i});
        end
        
    end
    
    for i = 1:numel(fnames)
        
        if strncmpi(fliplr(fnames{i}), 'loV', 3)
            design = rmiffield(design, fnames{i});
        end
        
    end
    
end
function [FemmProblem, outermagsep, coillabellocs] = slottedLfemmprob_torus(design, varargin)
% creates a FemmProblem structure for a slotted torus axial flux permanent
% magnet machine

    magsep = 2*design.g + 2*design.tc + 2*design.tsb + design.ty;
    
    Inputs.NStages = 1;
    Inputs.MagArrangement = 'NN';
    Inputs.NWindingLayers = 1;
     % set a suitible current for the inductance simulation in the circuit
    % for phase 1
    Inputs.CoilCurrent = inductancesimcurrent(design.CoilArea, design.CoilTurns);
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.Rmo - design.Rmi);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = 0;
    Inputs.MagnetSpaceGroup = 0;
    Inputs.MagnetSpaceMaterial = 1;
    Inputs.CoilGroup = 0;
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm(design.tm, design.taumm, 1/40);
    
    if min(design.tbi > 0)
        Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(min(design.tbi), design.taupm, 1/40);
    else
        Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(max(design.tbi), design.taupm, 1/40);
    end
    
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, design.taumm, 1/10), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm(magsep, design.taupm, 1/50);
    Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(design.tsg, design.tausgm, 1/50);
    Inputs.YokeRegionMeshSize = min( choosemesharea_mfemm(design.ty, 2*design.taupm, 1/40), ...
                                     choosemesharea_mfemm(design.tc, design.taucs, 1/40)  );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm(design.tc, design.taucs);
    Inputs.Tol = 1e-5;
    Inputs.NSlots = 2*design.Phases;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    slotsperpole = design.Qs / design.Poles;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    [FemmProblem, matinds] = addmaterials_mfemm(FemmProblem, ...
        {design.MagFEASimMaterials.Magnet, design.MagFEASimMaterials.FieldBackIron, design.MagFEASimMaterials.ArmatureYoke, design.MagFEASimMaterials.ArmatureCoil});
                 
%     MagnetMatInd = matinds(1);
    BackIronMatInd = matinds(2);
    YokeMatInd = matinds(3);
    CoilMatInd = matinds(4);
    
    if Inputs.MagnetSpaceMaterial ~= 1
        warning('SLOTTEDLFEMMPROB_TORUS:magspacenotair', ...
            ['The material in the region between the magnets is not air, ', ...
             'this may result in errors in the results due to the way the ', ...
             'inductance simulation is drawn']);
    end
    
    % draw the torus rotor according to the spec in the design strucure
    [FemmProblem, outermagsep, innerstagewidth] = torusrotor2dfemmprob( ...
        design.taupm, design.taumm, design.tm, design.tbi, magsep, ...
        'NStages', Inputs.NStages, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', 1, ... % don't add magnet label
        'BackIronMaterial', BackIronMatInd, ...
        'OuterRegionsMaterial', 1, ... % Air
        'MagnetSpaceMaterial', Inputs.MagnetSpaceMaterial, ... % don't add magnet space labels
        'MagnetGroup', Inputs.MagnetGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'Tol', Inputs.Tol);
    
    newupperbound = Inputs.NSlots*design.taupm/slotsperpole;
    scalefactor = newupperbound / (2*design.taupm);
    
    % find all nodes that are not at the very top and bottom and remove
    % them. The segments they are attached to will be removed as well
%     nodeids = [];
    for i = 1:numel(FemmProblem.Nodes)
%         if FemmProblem.Nodes(i).Coords(2) > 0 && ...
%            FemmProblem.Nodes(i).Coords(2) < (2*design.taupm - eps(2*design.taupm))
%             
%             nodeids = [nodeids, i - 1];
%             
%         end

        FemmProblem.Nodes(i).Coords(2) = FemmProblem.Nodes(i).Coords(2) * scalefactor;
    end
                
%     FemmProblem = removenodes_mfemm(FemmProblem, nodeids);
%     
%     % move all remaining nodes at a position greater than zero so that they
%     % are at the top of the slots
%     newypos = Inputs.NSlots*design.taupm/slotsperpole;
%     
%     for i = 1:numel(FemmProblem.Nodes)
%         if FemmProblem.Nodes(i).Coords(2) > 0
%             FemmProblem.Nodes(i).Coords(2) = newypos;
%         end
%     end
%     
%     newlabypos = newypos / 2;
    
    for i = 1:numel(FemmProblem.BlockLabels)
%         if FemmProblem.BlockLabels(i).Coords(2) > 0
%             FemmProblem.BlockLabels(i).Coords(2) = newlabypos;
%         end

        FemmProblem.BlockLabels(i).Coords(2) = FemmProblem.BlockLabels(i).Coords(2) * scalefactor;
        
    end
    
    
    % draw the stator slots for all stages
    [FemmProblem, yokenodeids, coillabellocs] = axialfluxinnerstator2dfemmprob( ...
        innerstagewidth, design.Qs, design.Poles, design.taupm, design.taucs, ...
        design.tausgm, design.ty, design.tc, design.tsb, design.tsg, ...
        'NStators', Inputs.NStages, ...
        'NWindingLayers', Inputs.NWindingLayers, ...
        'FemmProblem', FemmProblem, ...
        'ToothMaterial', YokeMatInd, ...
        'ToothRegionMeshSize', Inputs.YokeRegionMeshSize, ...
        'ShoeGapMaterial', 1, ...
        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
        'Tol', Inputs.Tol, ...
        'NSlots', Inputs.NSlots);
    
    % link the rotor stages along the top and bottom, add antiperiodic
    % boundaries, and add the coil regions.

    gapedgenodes = [-outermagsep/2, 0;
                    -outermagsep/2+magsep, 0;
                    -outermagsep/2+magsep, newupperbound;
                    -outermagsep/2, newupperbound];
                
    FemmProblem = slottedcommonfemmprob_torus(FemmProblem, design, ...
        Inputs, magsep, gapedgenodes, innerstagewidth, coillabellocs, ...
        yokenodeids, newupperbound/2, BackIronMatInd, YokeMatInd, CoilMatInd);                

end

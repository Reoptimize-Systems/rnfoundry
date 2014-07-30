function [FemmProblem, coillabellocs, inslabellocs] = slottedLfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotted radial flux permanent
% magnet machine for an inductance simulation
%
% Syntax
%
% [FemmProblem, coillabellocs] = slottedLfemmprob_radial(design, varargin)
%
% 

    design.thetas = 2*pi/design.Qs;
    
    Inputs.ArmatureType = 'external';
    Inputs.MagArrangement = 'NN';
    Inputs.NWindingLayers = 1;
     % set a suitible current for the inductance simulation in the circuit
    % for phase 1
    Inputs.CoilCurrent = inductancesimcurrent(annularsecarea(design.Rci, design.Rco, design.thetac), ...
                                              design.CoilTurns);
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.ls);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = 0;
    Inputs.MagnetSpaceGroup = 0;
    Inputs.RotorBackIronGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.ShoeGroup = 0;
    Inputs.ArmatureBackIronGroup = 0;
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm(design.tm, design.Rmm*design.thetam, 1/40);
    
    if min(design.tbi > 0)
        Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(min(design.tbi), 2*design.Rbm*design.thetap, 1/40);
    else
        Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(max(design.tbi), 2*design.Rbm*design.thetap, 1/40);
    end
    
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, design.Rmm*design.thetam, 1/10), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm(design.g, design.Rmm*design.thetap, 1/50);
    Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(design.tsg, design.Rcm*(design.thetas-mean(design.thetac))/2, 1/50);
    Inputs.YokeRegionMeshSize = min( choosemesharea_mfemm(design.ty, 2*design.Rym*design.thetap, 1/40), ...
                                     choosemesharea_mfemm(design.tc, design.Rcm*mean(design.thetac), 1/40)  );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm(design.tc, design.Rcm*mean(design.thetac));
    Inputs.Tol = 1e-5;
    Inputs.NSlots = 2*design.Phases;
    Inputs.DrawOuterRegions = true;
    Inputs.CoilInsRegionMeshSize = -1;
    Inputs.DrawCoilInsulation = false;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    slotsperpole = design.Qs / design.Poles;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    [FemmProblem, matinds] = addmaterials_mfemm(FemmProblem, ...
        {design.MagFEASimMaterials.AirGap, design.MagFEASimMaterials.Magnet, design.MagFEASimMaterials.FieldBackIron, ...
         design.MagFEASimMaterials.ArmatureYoke, design.MagFEASimMaterials.ArmatureCoil});
    
    GapMatInd = matinds(1);
%     MagnetMatInd = matinds(1);
    BackIronMatInd = matinds(2);
    YokeMatInd = matinds(3);
    CoilMatInd = matinds(4);
    
    if GapMatInd ~= 1
        warning('SLOTTEDLFEMMPROB_TORUS:magspacenotair', ...
            ['The material in the region between the magnets is not air, ', ...
             'this may result in errors in the results due to the way the ', ...
             'inductance simulation is drawn']);
    end
    
    switch Inputs.ArmatureType
        case 'external'
            % single inner facing stator
            drawnrotors = [false, true];
            rrotor = design.Rmo;
            drawnstatorsides = [1, 0];
            Rs = design.Rmo + design.g + design.tc + design.tsb + design.ty/2;
        case 'internal'
            % single outer facing stator
            drawnrotors = [true, false];
            rrotor = design.Rmo;
            drawnstatorsides = [0, 1]; 
            Rs = design.Rmi - design.g - design.tc - design.tsb - design.ty/2;
        case 'di'
            % double internal stator (mags on outside)
%             drawnrotors = [true, true];
%             rrotor = [ design.Rmo, design.Rmo + 2* (design.g + design.tc + design.ty/2) ];
%             drawnstatorsides = [1, 1];
%             Rs = design.Rmo(1) + design.g + design.tc + design.ty/2;
            error('not yet supported');
        case 'do'
            % double outer/external stator (mags on inside)
            error('not yet supported');
            
        otherwise
            error('Unrecognised ArmatureType option.')
                
    end
    
    newupperbound = Inputs.NSlots*design.thetap/slotsperpole;
    scalefactor = newupperbound / (2*design.thetap);
    
    % create a modified copy of the design 
    Ldesign = design;
    Ldesign.thetap = Ldesign.thetap * scalefactor;
    Ldesign.thetam = Ldesign.thetam * scalefactor;
    
    % draw the rotor according to the spec in the design strucure
    [FemmProblem, magcornerids, linktb] = radialfluxrotor2dfemmprob( ...
        Ldesign.thetap, Ldesign.thetam, design.tm, design.tbi, ...
        drawnrotors, rrotor, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', 1, ... % don't add magnet material for inductance, always air
        'BackIronMaterial', BackIronMatInd, ...
        'OuterRegionsMaterial', GapMatInd, ...
        'MagnetSpaceMaterial', GapMatInd, ... % don't add magnet space labels
        'MagnetGroup', Inputs.MagnetGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'BackIronGroup', Inputs.RotorBackIronGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'Tol', Inputs.Tol, ...
        'DrawOuterRegions', Inputs.DrawOuterRegions );


    if numel (design.tc) > 1
        coilbasefrac = design.tc(2) / design.tc(1);
    else
        coilbasefrac = 0.05;
    end
    
    if isfield (design, 'ShoeCurveControlFrac')
        shoecurvefrac = design.ShoeCurveControlFrac;
    else
        shoecurvefrac = 0.5;
    end
    
    % draw the stator slots for all stages
    [FemmProblem, yokenodeids, coillabellocs, inslabellocs] = radialfluxstator2dfemmprob ( ...
        design.Qs, design.Poles, Rs, design.thetap, design.thetac, ...
        design.thetasg, design.ty, design.tc, design.tsb, design.tsg, drawnstatorsides, ...
        'NWindingLayers', Inputs.NWindingLayers, ...
        'FemmProblem', FemmProblem, ...
        'ShoeGapMaterial', GapMatInd, ...
        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
        'ShoeGroup', Inputs.ArmatureBackIronGroup, ...
        'Tol', Inputs.Tol, ...
        'NSlots', Inputs.NSlots, ...
        'DrawCoilInsulation', Inputs.DrawCoilInsulation, ...
        'CoilInsulationThickness', design.CoilInsulationThickness, ...
        'CoilBaseFraction', coilbasefrac, ...
        'ShoeCurveControlFrac', shoecurvefrac );
    
    
    % link the rotor stages along the top and bottom, add antiperiodic
    % boundaries, and add the coil regions.
    FemmProblem = slottedcommonfemmprob_radial ( FemmProblem, ...
                                                 Ldesign, ...
                                                 Inputs, ...
                                                 magcornerids, ...
                                                 Rs, ...
                                                 coillabellocs, ...
                                                 inslabellocs, ...
                                                 yokenodeids, ...
                                                 Ldesign.thetap, ...
                                                 BackIronMatInd, ...
                                                 YokeMatInd, ...
                                                 CoilMatInd, ...
                                                 GapMatInd, ...
                                                 linktb );             

end

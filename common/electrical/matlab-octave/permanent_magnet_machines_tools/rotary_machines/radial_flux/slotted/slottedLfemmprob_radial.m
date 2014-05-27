function [FemmProblem, coillabellocs] = slottedLfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotted torus axial flux permanent
% magnet machine

    design.thetas = 2*pi/design.Qs;
    
    Inputs.StatorType = 'si';
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
    Inputs.MagnetSpaceMaterial = 1;
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
    Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(design.tsg, design.Rcm*(design.thetas-design.thetac)/2, 1/50);
    Inputs.YokeRegionMeshSize = min( choosemesharea_mfemm(design.ty, 2*design.Rym*design.thetap, 1/40), ...
                                     choosemesharea_mfemm(design.tc, design.Rcm*design.thetac, 1/40)  );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm(design.tc, design.Rcm*design.thetac);
    Inputs.Tol = 1e-5;
    Inputs.NSlots = 2*design.Phases;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    slotsperpole = design.Qs / design.Poles;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    [FemmProblem, matinds] = addmaterials_mfemm(FemmProblem, ...
        {design.MagnetMaterial, design.BackIronMaterial, design.YokeMaterial, design.CoilMaterial});
                 
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
    
    switch Inputs.StatorType
        case 'si'
            % single inner facing stator
            drawnrotors = [false, true];
            rrotor = design.Rmo;
            drawnstatorsides = [1, 0];
            Rs = design.Rmo + design.g + design.tc + design.tsb + design.ty/2;
        case 'so'
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
            error('Unrecognised StatorType option.')
                
    end
    
    newupperbound = Inputs.NSlots*design.thetap/slotsperpole;
    scalefactor = newupperbound / (2*design.thetap);
    
    % create a modified copy of the design 
    Ldesign = design;
    Ldesign.thetap = Ldesign.thetap * scalefactor;
    Ldesign.thetam = Ldesign.thetam * scalefactor;
    
    % draw the torus rotor according to the spec in the design strucure
    [FemmProblem, magcornerids, linktb] = radialfluxrotor2dfemmprob( ...
        Ldesign.thetap, Ldesign.thetam, design.tm, design.tbi, ...
        drawnrotors, rrotor, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', 1, ... % don't add magnet label
        'BackIronMaterial', BackIronMatInd, ...
        'OuterRegionsMaterial', 1, ... % Air
        'MagnetSpaceMaterial', Inputs.MagnetSpaceMaterial, ... % don't add magnet space labels
        'MagnetGroup', Inputs.MagnetGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'BackIronGroup', Inputs.RotorBackIronGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'Tol', Inputs.Tol);
    
    % find all nodes that are not at the very top and bottom and remove
    % them. The segments they are attached to will be removed as well
% %     nodeids = [];
%     for i = 1:numel(FemmProblem.Nodes)
% %         if FemmProblem.Nodes(i).Coords(2) > 0 && ...
% %            FemmProblem.Nodes(i).Coords(2) < (2*design.taupm - eps(2*design.taupm))
% %             
% %             nodeids = [nodeids, i - 1];
% %             
% %         end
% 
%         FemmProblem.Nodes(i).Coords(2) = FemmProblem.Nodes(i).Coords(2) * scalefactor;
%     end
                
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
    
%     for i = 1:numel(FemmProblem.BlockLabels)
% %         if FemmProblem.BlockLabels(i).Coords(2) > 0
% %             FemmProblem.BlockLabels(i).Coords(2) = newlabypos;
% %         end
% 
%         FemmProblem.BlockLabels(i).Coords(2) = FemmProblem.BlockLabels(i).Coords(2) * scalefactor;
%         
%     end
    
    
    % draw the stator slots for all stages
    [FemmProblem, yokenodeids, coillabellocs] = radialfluxstator2dfemmprob( ...
        design.Qs, design.Poles, Rs, design.thetap, design.thetac, ...
        design.thetasg, design.ty, design.tc, design.tsb, design.tsg, drawnstatorsides, ...
        'NWindingLayers', Inputs.NWindingLayers, ...
        'FemmProblem', FemmProblem, ...
        ...'ToothMaterial', YokeMatInd, ...
        ...'ToothRegionMeshSize', Inputs.YokeRegionMeshSize, ...
        'ShoeGapMaterial', 1, ...
        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
        'ShoeGroup', Inputs.ArmatureBackIronGroup, ...
        'Tol', Inputs.Tol, ...
        'NSlots', Inputs.NSlots);
    
    
    % link the rotor stages along the top and bottom, add antiperiodic
    % boundaries, and add the coil regions.
    FemmProblem = slottedcommonfemmprob_radial( FemmProblem, ...
                                                Ldesign, ...
                                                Inputs, ...
                                                magcornerids, ...
                                                Rs, ...
                                                coillabellocs, ...
                                                yokenodeids, ...
                                                Ldesign.thetap, ...
                                                BackIronMatInd, ...
                                                YokeMatInd, ...
                                                CoilMatInd, ...
                                                Inputs.ArmatureBackIronGroup, ...
                                                linktb );             

end

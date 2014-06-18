function [FemmProblem, coillabellocs, yokenodeids] = slottedfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotted radial flux permanent
% magnet machine
%
% Syntax
%
% [FemmProblem, coillabellocs, yokenodeids] = slottedfemmprob_radial(design)
% [FemmProblem, coillabellocs, yokenodeids] = slottedfemmprob_radial(..., 'Parameter', Value)
%

    % First set up some default inputs
    Inputs.ArmatureType = 'external';
    Inputs.NWindingLayers = 1;
    Inputs.CoilCurrent = 0;
    Inputs.MagArrangement = 'NN';
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.ls, 'MinAngle', 15);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup = [];
    Inputs.RotorBackIronGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.ShoeGroup = 0;
    Inputs.ArmatureBackIronGroup = 0;
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm(design.tm, (design.Rmm*design.thetam), 1/10);
    Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(min(design.tbi), 2*(design.Rbm*design.thetap), 1/10);
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm(design.g, (design.Rmm*design.thetap), 1/10);
    Inputs.DrawOuterRegions = true;
    
    if design.tsg > 1e-5
        if design.tsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20);
        end
    else
        if design.tsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = -1;
        end
    end
    Inputs.YokeRegionMeshSize = mean( [choosemesharea_mfemm(design.ty, 2*(design.Rym*design.thetap), 1/10), ...
                                       choosemesharea_mfemm(design.tc, (design.Rcm*(design.thetas-design.thetac)), 1/10)] );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm(design.tc, (design.Rcm*design.thetac));
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(design.thetap, ...
                                     Inputs.Position, ...
                                     Inputs.FractionalPolePosition, ...
                                     Inputs.RotorAnglePosition);
    
    Inputs.NSlots = 2*design.Qs/design.Poles;
    
    FemmProblem = Inputs.FemmProblem;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    [FemmProblem, matinds] = addmaterials_mfemm (FemmProblem, ...
        {design.GapMaterial, design.MagnetMaterial, design.BackIronMaterial, ...
         design.YokeMaterial, design.CoilMaterial} );
                 
    GapMatInd = matinds(1);
    MagnetMatInd = matinds(2);
    BackIronMatInd = matinds(3);
    YokeMatInd = matinds(4);
    CoilMatInd = matinds(5);
    
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
            rrotor = design.Rmi;
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

    % draw the radial rotor according to the spec in the design strucure
    [FemmProblem, magcornerids, linktb] = radialfluxrotor2dfemmprob( ...
        design.thetap, design.thetam, design.tm, design.tbi, drawnrotors, rrotor, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', MagnetMatInd, ...
        'BackIronMaterial', BackIronMatInd, ...
        'OuterRegionsMaterial', GapMatInd, ... % ususally Air
        'MagnetSpaceMaterial', GapMatInd, ... % usually Air
        'MagnetGroup', Inputs.MagnetGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'BackIronGroup', Inputs.RotorBackIronGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'Tol', Inputs.Tol);
    
    % draw the stator slots
    [FemmProblem, yokenodeids, coillabellocs] = radialfluxstator2dfemmprob( ...
        design.Qs, design.Poles, Rs, design.thetap, design.thetac, ...
        design.thetasg, design.ty, design.tc, design.tsb, design.tsg, drawnstatorsides, ...
        'NWindingLayers', Inputs.NWindingLayers, ...
        'FemmProblem', FemmProblem, ...
        ... 'ToothMaterial', YokeMatInd, ...
        ... 'ToothRegionMeshSize', Inputs.YokeRegionMeshSize, ...
        'ShoeGapMaterial', GapMatInd, ...
        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
        'ShoeGroup', Inputs.ArmatureBackIronGroup, ...
        'Tol', Inputs.Tol);

    % Complete the design using the common radial drawing function
    FemmProblem = slottedcommonfemmprob_radial( FemmProblem, ...
                                                design, ...
                                                Inputs, ...
                                                magcornerids, ...
                                                Rs, ...
                                                coillabellocs, ...
                                                yokenodeids, ...
                                                design.thetap, ...
                                                BackIronMatInd, ...
                                                YokeMatInd, ...
                                                CoilMatInd, ...
                                                Inputs.ArmatureBackIronGroup, ...
                                                linktb, ...
                                                GapMatInd, ...
                                                Inputs.DrawOuterRegions );

end

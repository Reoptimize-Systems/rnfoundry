function [FemmProblem, coillabellocs] = slottedfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotted radial flux permanent
% magnet machine
%
% Syntax
%
% [FemmProblem, coillabellocs] = slottedfemmprob_radial(design)
% [FemmProblem, coillabellocs] = slottedfemmprob_radial(..., 'Parameter', Value)
%

    % First set up some default inputs
    Inputs.ArmatureType = 'external';
    Inputs.NWindingLayers = 1;
    Inputs.CoilCurrent = 0;
    Inputs.MagArrangement = 'NN';
    Inputs.PolarisationType = 'constant';
    if isfield (design, 'MagnetPolarisation') && ischar (design.MagnetPolarisation)
        Inputs.PolarisationType = design.MagnetPolarisation;
    end
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.ls, 'MinAngle', 15);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup = [];
    Inputs.RotorBackIronGroup = [];
    Inputs.RotorOuterRegionGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.ArmatureBackIronGroup = [];
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm(design.tm, (design.Rmm*design.thetam), 1/10);
    Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(min(design.tbi), 2*(design.Rbm*design.thetap), 1/10);
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm(design.g, (design.Rmm*design.thetap), 1/10);
    Inputs.DrawOuterRegions = true;
    Inputs.DrawCoilInsulation = false;
    Inputs.CoilInsRegionMeshSize = -1;
    
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
                                       choosemesharea_mfemm(design.tc(1), (design.Rcm*(design.thetas-mean(design.thetac))), 1/10)] );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm(design.tc(1), (design.Rcm*mean(design.thetac)));
    Inputs.Tol = 1e-5;
    Inputs.SimType = 'Magnetics';
    Inputs.MaterialsLibrary = '';
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    if isempty(Inputs.ArmatureBackIronGroup)
        [FemmProblem, Inputs.ArmatureBackIronGroup] = addgroup_mfemm(FemmProblem, 'ArmatureBackIron');
    end
    
    if isempty (Inputs.MaterialsLibrary)
        if strncmpi (Inputs.SimType, 'Magnetics', 1)
            FemmProblem.ProbInfo.Domain = 'Magnetics';
            Inputs.MaterialsLibrary = fullfile(fileparts (which ('matstr2matstruct_mfemm')), '..', 'matlib.mat');
        elseif strncmpi (Inputs.SimType, 'HeatFlow', 1)
            FemmProblem.ProbInfo.Domain = 'HeatFlow';
            Inputs.MaterialsLibrary = fullfile(fileparts (which ('matstr2matstruct_mfemm')), '..', 'heatlib.mat');
        else
            error ('Unrecognised SimType');
        end
    end
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(design.thetap, ...
                                     Inputs.Position, ...
                                     Inputs.FractionalPolePosition, ...
                                     Inputs.RotorAnglePosition);
    
    Inputs.NSlots = 2*design.Qs/design.Poles;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    if strncmpi (Inputs.SimType, 'Magnetics', 1)
    
        [FemmProblem, matinds] = addmaterials_mfemm (FemmProblem, ...
            { design.MagFEASimMaterials.AirGap, ...
              design.MagFEASimMaterials.Magnet, ...
              design.MagFEASimMaterials.FieldBackIron, ...
              design.MagFEASimMaterials.ArmatureYoke, ...
              design.MagFEASimMaterials.ArmatureCoil }, ...
             'MaterialsLibrary', Inputs.MaterialsLibrary );
         
    elseif strncmpi (Inputs.SimType, 'HeatFlow', 1)
    
        [FemmProblem, matinds] = addmaterials_mfemm (FemmProblem, ...
            { design.HeatFEASimMaterials.AirGap, ...
              design.HeatFEASimMaterials.Magnet, ...
              design.HeatFEASimMaterials.FieldBackIron, ...
              design.HeatFEASimMaterials.ArmatureYoke, ...
              design.HeatFEASimMaterials.ArmatureCoil }, ...
             'MaterialsLibrary', Inputs.MaterialsLibrary );
    
    end
                 
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
            Rs = design.Rmo + design.g + design.tc(1) + design.tsb + design.ty/2;
        case 'internal'
            % single outer facing stator
            drawnrotors = [true, false];
            rrotor = design.Rmi;
            drawnstatorsides = [0, 1]; 
            Rs = design.Rmi - design.g - design.tc(1) - design.tsb - design.ty/2;
        case 'di'
            % double internal stator (mags on outside)
%             drawnrotors = [true, true];
%             rrotor = [ design.Rmo, design.Rmo + 2* (design.g + design.tc(1) + design.ty/2) ];
%             drawnstatorsides = [1, 1];
%             Rs = design.Rmo(1) + design.g + design.tc(1) + design.ty/2;
            error('not yet supported');
        case 'do'
            % double outer/external stator (mags on inside)
            error('not yet supported');
            
        otherwise
            error('Unrecognised ArmatureType option.')
                
    end

    % draw the radial rotor according to the spec in the design strucure
    [FemmProblem, magcornerids, linktb] = radialfluxrotor2dfemmprob ( ...
        design.thetap, design.thetam, design.tm, design.tbi, drawnrotors, rrotor, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'PolarisationType', Inputs.PolarisationType, ...
        'MagnetMaterial', MagnetMatInd, ...
        'BackIronMaterial', BackIronMatInd, ...
        'OuterRegionsMaterial', GapMatInd, ... % ususally Air
        'MagnetSpaceMaterial', GapMatInd, ... % usually Air
        'MagnetGroup', Inputs.MagnetGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'BackIronGroup', Inputs.RotorBackIronGroup, ...
        'OuterRegionGroup', Inputs.RotorOuterRegionGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'DrawOuterRegions', Inputs.DrawOuterRegions, ...
        'Tol', Inputs.Tol );
        
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
    
    % draw the stator slots
    [FemmProblem, yokenodeids, coillabellocs, inslabellocs] = radialfluxstator2dfemmprob ( ...
        design.Qs, design.Poles, Rs, design.thetap, design.thetac, ...
        design.thetasg, design.ty, design.tc(1), design.tsb, design.tsg, drawnstatorsides, ...
        'NWindingLayers', Inputs.NWindingLayers, ...
        'FemmProblem', FemmProblem, ...
        'ShoeGapMaterial', GapMatInd, ...
        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
        'ArmatureIronGroup', Inputs.ArmatureBackIronGroup, ...
        'Tol', Inputs.Tol, ...
        'DrawCoilInsulation', Inputs.DrawCoilInsulation, ...
        'CoilInsulationThickness', design.CoilInsulationThickness, ...
        'CoilBaseFraction', coilbasefrac, ...
        'ShoeCurveControlFrac', shoecurvefrac );

    % Complete the design using the common radial drawing function
    FemmProblem = slottedcommonfemmprob_radial( FemmProblem, ...
                                                design, ...
                                                Inputs, ...
                                                magcornerids, ...
                                                Rs, ...
                                                coillabellocs, ...
                                                inslabellocs, ...
                                                yokenodeids, ...
                                                design.thetap, ...
                                                BackIronMatInd, ...
                                                YokeMatInd, ...
                                                CoilMatInd, ...
                                                GapMatInd, ...
                                                linktb );
                                                

end

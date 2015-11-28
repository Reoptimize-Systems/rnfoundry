function [FemmProblem, rotorinfo, statorinfo] = slottedfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotted radial flux permanent
% magnet machine
%
% Syntax
%
% [FemmProblem, coillabellocs] = slottedfemmprob_radial (design)
% [FemmProblem, coillabellocs] = slottedfemmprob_radial (..., 'Parameter', Value)
%
% 
% Inputs
%
%  design - Structure containing the design specification.
%
%  Finer control over the drawing can be made using a number of optional
%  arguments supplied as parameter-value pairs, i.e. a string followed by
%  the desired value. The optional parameters are described below:
%
%
%   'DrawingType' - string describing the desired type of simulation
%       drawing. Currently there is only one option '2PoleMagnetRotation', also
%       the default. With this option the magnets are drwan rotated to the
%       desired position. See the 'Position' parameter below for more
%       details.
%
%   'Position' = 0;
%
%   'NBoundaryPositions' = 10;
%
%   'ArmatureType' = 'external'
%
%   'NWindingLayers' = 1;
%
%   'CoilCurrent' = zeros (1,design.Phases)
%
%   'MagArrangement' = 'NN'
%
%   'PolarisationType' = 'constant'
%
%   'PolarisationType' = design.MagnetPolarisation
%
%   'FemmProblem' = newproblem_mfemm('planar', 'Depth', design.ls, 'MinAngle', 15)
%
%   'FractionalPolePosition' = [];
%
%   'RotorAnglePosition' = [];
%
%   'MagnetGroup' = [];
%
%   'MagnetSpaceGroup' = [];
%
%   'RotorBackIronGroup' = [];
%
%   'RotorOuterRegionGroup' = [];
%
%   'CoilGroup' = 0;
%
%   'ArmatureBackIronGroup' = [];
%
%   'MagnetRegionMeshSize' = choosemesharea_mfemm(design.tm, (design.Rmm*design.thetam), 1/10);
%
%   'BackIronRegionMeshSize' = choosemesharea_mfemm(min(design.tbi), 2*(design.Rbm*design.thetap), 1/10);
%
%   'OuterRegionsMeshSize' = [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1];
%
%   'AirGapMeshSize' = choosemesharea_mfemm(design.g, (design.Rmm*design.thetap), 1/10);
%
%   'DrawOuterRegions' = true;
%
%   'DrawCoilInsulation' = false;
%
%   'CoilInsRegionMeshSize' = -1;
%
%   'ShoeGapRegionMeshSize' = choosemesharea_mfemm(max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20);
%
%   'YokeRegionMeshSize'
%
%   'CoilRegionMeshSize'
%
%   'Tol' = 1e-5;
%
%   'SimType' = 'Magnetics';
%
%   'MaterialsLibrary' = '';
%
%   'NPolePairs' = 1;
%
%   'NSlots' = NPolePairs*2*design.Qs/design.Poles;
%

    % First set up some default inputs
    Inputs.DrawingType = 'MagnetRotation';
    Inputs.NBoundaryPositions = 10;
    Inputs.BoundaryShift = 0;
    Inputs.ArmatureType = 'external';
    Inputs.NWindingLayers = nan;
    Inputs.CoilCurrent = zeros (1,design.Phases);
    Inputs.MagArrangement = 'NN';
    Inputs.PolarisationType = 'constant';
    if isfield (design, 'MagnetPolarisation') && ischar (design.MagnetPolarisation)
        Inputs.PolarisationType = design.MagnetPolarisation;
    end
    Inputs.FemmProblem = newproblem_mfemm ('planar', 'Depth', design.ls, 'MinAngle', 15);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup = [];
    Inputs.RotorBackIronGroup = [];
    Inputs.RotorOuterRegionGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.ArmatureBackIronGroup = [];
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm (design.tm, (design.Rmm*design.thetam), 1/10);
    Inputs.BackIronRegionMeshSize = choosemesharea_mfemm (min(design.tbi), 2*(design.Rbm*design.thetap), 1/10);
    Inputs.RotorOuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1];
    Inputs.StatorOuterRegionsMeshSize = [];
    Inputs.StatorOuterRegionMaterials = {};
    Inputs.AirGapMeshSize = choosemesharea_mfemm (design.g, (design.Rmm*design.thetap), 1/10);
    Inputs.DrawOuterRegions = true;
    Inputs.StatorOuterRegionSize = [];
    Inputs.DrawCoilInsulation = false;
    Inputs.CoilInsRegionMeshSize = -1;
    
    if design.tsg > 1e-5
        if design.tsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (design.tsb, (design.Rmo*design.thetasg), 1/20);
        end
    else
        if design.tsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (design.tsb, (design.Rmo*design.thetasg), 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = -1;
        end
    end
    Inputs.YokeRegionMeshSize = mean( [choosemesharea_mfemm(design.ty, 2*(design.Rym*design.thetap), 1/10), ...
                                       choosemesharea_mfemm(design.tc(1), (design.Rcm*(design.thetas-mean(design.thetac))), 1/10)] );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm (design.tc(1), (design.Rcm*mean(design.thetac)));
    Inputs.Tol = 1e-5;
    Inputs.SimType = 'Magnetics';
    Inputs.MaterialsLibrary = '';
    Inputs.NPolePairs = 1;
    
    Inputs = parseoptions (Inputs, varargin);
    
    Inputs.NSlots = Inputs.NPolePairs*2*design.Qs/design.Poles;
    
    if isnan(Inputs.NWindingLayers)
        if isfield (design, 'CoilLayers')
            Inputs.NWindingLayers = design.CoilLayers;
        else
            Inputs.NWindingLayers = 1;
            warning ('Number of winding layers not specified, using 1.');
        end
    end
    
    FemmProblem = Inputs.FemmProblem;
    
    if isempty (Inputs.ArmatureBackIronGroup) ...
            && ~isfield (FemmProblem.Groups, 'ArmatureBackIron')
        [FemmProblem, Inputs.ArmatureBackIronGroup] = addgroup_mfemm (FemmProblem, 'ArmatureBackIron');
    end
    
    if isempty (Inputs.MaterialsLibrary)
        if strncmpi (Inputs.SimType, 'Magnetics', 1)
            FemmProblem.ProbInfo.Domain = 'Magnetics';
            Inputs.MaterialsLibrary = fullfile (fileparts (which ('matstr2matstruct_mfemm')), '..', 'matlib.mat');
        elseif strncmpi (Inputs.SimType, 'HeatFlow', 1)
            FemmProblem.ProbInfo.Domain = 'HeatFlow';
            Inputs.MaterialsLibrary = fullfile (fileparts (which ('matstr2matstruct_mfemm')), '..', 'heatlib.mat');
        else
            error ('Unrecognised SimType');
        end
    end
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos ( design.thetap, ...
                                       Inputs.Position, ...
                                       Inputs.FractionalPolePosition, ...
                                       Inputs.RotorAnglePosition );
    
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
            outerR = design.Ryo;
        case 'internal'
            % single outer facing stator
            drawnrotors = [true, false];
            rrotor = design.Rmi;
            drawnstatorsides = [0, 1]; 
            Rs = design.Rmi - design.g - design.tc(1) - design.tsb - design.ty/2;
            outerR = design.Rbo;
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
    
    lastslot = 1;
    XShift = 0;
    YShift = 0;
%     coillabellocs = [];
    tbboundseginds = [];
    
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
    
    % define the block properties of the core region
    yokeBlockProps.BlockType = FemmProblem.Materials(YokeMatInd).Name;
    yokeBlockProps.MaxArea = Inputs.BackIronRegionMeshSize;
    yokeBlockProps.InCircuit = '';
    yokeBlockProps.InGroup = Inputs.ArmatureBackIronGroup;
    
    switch Inputs.DrawingType
        
        case 'MagnetRotation'
    
            % draw the radial rotor according to the spec in the design strucure
            [FemmProblem, rotorinfo] = radialfluxrotor2dfemmprob ( ...
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
                'OuterRegionsMeshSize', Inputs.RotorOuterRegionsMeshSize, ...
                'Position', Inputs.Position, ...
                'DrawOuterRegions', Inputs.DrawOuterRegions, ...
                'NPolePairs', Inputs.NPolePairs, ...
                'Tol', Inputs.Tol, ...
                'YShift', YShift );


            % draw the stator slots
            [FemmProblem, statorinfo] = radialfluxstator2dfemmprob ( ...
                design.Qs, design.Poles, Rs, design.thetap, design.thetac, ...
                design.thetasg, design.ty, design.tc(1), design.tsb, design.tsg, drawnstatorsides, ...
                'NWindingLayers', Inputs.NWindingLayers, ...
                'FemmProblem', FemmProblem, ...
                'ShoeGapMaterial', GapMatInd, ...
                'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                ...'ArmatureIronGroup', Inputs.ArmatureBackIronGroup, ...
                'Tol', Inputs.Tol, ...
                'DrawCoilInsulation', Inputs.DrawCoilInsulation, ...
                'CoilInsulationThickness', design.CoilInsulationThickness, ...
                'CoilBaseFraction', coilbasefrac, ...
                'ShoeCurveControlFrac', shoecurvefrac, ...
                'NSlots', Inputs.NSlots, ...
                'YShift', YShift );

%             coillabellocs = [coillabellocs; statorinfo.CoilLabelLocations];

            % Complete the design using the common radial drawing function
            Inputs.AddAllCoilsToCircuits = true;
%             Inputs.StartSlot = lastslot;
            
            [FemmProblem, statorinfo] = stator_iron_boundary (FemmProblem, design, Inputs, statorinfo, Rs, XShift, YShift);
            
            [FemmProblem, statorinfo] = stator_outer_regions (FemmProblem, design, Inputs, statorinfo, Rs, GapMatInd, XShift, YShift);
            
            if abs(tau-(design.thetap * 2 * Inputs.NPolePairs)) > Inputs.Tol
                % create boundaries for the air gap
                
                % add segments with periodic boundaries on the outer parts
                [FemmProblem, ~, airgapboundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Air Gap Periodic', 4);

                % bottom gap corner to bottom core corner
                [FemmProblem, commoninfo.BottomSegInds] = addsegments_mfemm (FemmProblem, ...
                                                   statorinfo.node_id_set_stator_inner(1), ...
                                                   rotorinfo.MagnetCornerIDs(1), ...
                                                   'BoundaryMarker', airgapboundname);
                
                % top gap corner to top core corner
                [FemmProblem, commoninfo.BottomSegInds(end+1)] = addsegments_mfemm (FemmProblem, ...
                                                   statorinfo.node_id_set_stator_inner(2), ...
                                                   rotorinfo.MagnetCornerIDs(2), ...
                                                   'BoundaryMarker', airgapboundname);
                
            end
            
%             [FemmProblem, commoninfo] = magnetrotation_complete_drawing (FemmProblem, design, Inputs, ...
%                                                            GapMatInd, ...
%                                                            rotorinfo.MagnetCornerIDs, ...
%                                                            Rs, ...
%                                                            statorinfo.OuterNodes, ...
%                                                            rotorinfo.LinkTopBottom, ...
%                                                            XShift, YShift );
                                                       
%             [FemmProblem, commoninfo] = slottedcommonfemmprob_radial ( FemmProblem, ...
%                                                         design, ...
%                                                         Inputs, ...
%                                                         rotorinfo.MagnetCornerIDs, ... magcornerids, ...
%                                                         Rs, ...
%                                                         statorinfo.CoilLabelLocations, ...
%                                                         statorinfo.InsulationLabelLocations, ...
%                                                         statorinfo.OuterNodes, ...
%                                                         design.thetap, ...
%                                                         BackIronMatInd, ...
%                                                         YokeMatInd, ...
%                                                         CoilMatInd, ...
%                                                         GapMatInd, ...
%                                                         rotorinfo.LinkTopBottom, ...
%                                                         XShift, ...
%                                                         YShift );

%             lastslot = lastslot + Inputs.NSlots;
%             YShift = YShift + (2.01 * outerR);
%             XShift = 0;
% 
%             tbboundseginds = [ [rotorinfo.TopSegInds, commoninfo.TopSegInds];
%                                [rotorinfo.BottomSegInds, commoninfo.BottomSegInds];
%                                tbboundseginds ];
% 
%             % rearrange the boundaries to link everything up correctly
% 
%             % first get all the boundary IDs in the appropriate segments
%             tbboundnames = getsegbounds (FemmProblem, tbboundseginds);
% 
%             % shift them round by one
%             tbboundnames = circshift (tbboundnames, [-1, 0]);
% 
%             % replace the boundaries
%             for ind = 1:numel (tbboundnames)
% 
%                 FemmProblem.Segments(tbboundseginds(ind)).BoundaryMarker = tbboundnames{ind};
% 
%             end
            
            % add the air gap labels
            switch Inputs.ArmatureType
        
                case 'external'

                    % Add block labels for the air gap
                    [labelloc(1),labelloc(2)]  = pol2cart (design.thetap, design.Rmo+design.g/2);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm ( labelloc(1,1), labelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize );

                    % add a block label for the yoke and teeth
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, Rs);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                                                         yokeBlockProps);

                case 'internal'

                    % Add block labels for the air gap
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, design.Rmi-design.g/2);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize);

                    % add a block label for the yoke and teeth
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, Rs);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                                                         yokeBlockProps);

                otherwise
                    error('Unrecognised ArmatureType option.')

            end
            
        case 'BoundarySwap'
            
            % draws the rotor and stator separately and links them at the
            % air gap using appropriate boundaries which can be swapped to
            % create relative motion without remeshing
            
            
            % draw the radial rotor according to the spec in the design strucure
            [FemmProblem, rotorinfo] = radialfluxrotor2dfemmprob ( ...
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
                'OuterRegionsMeshSize', Inputs.RotorOuterRegionsMeshSize, ...
                'Position', Inputs.Position, ...
                'DrawOuterRegions', Inputs.DrawOuterRegions, ...
                'NPolePairs', Inputs.NPolePairs, ...
                'Tol', Inputs.Tol, ...
                'YShift', 0 );
            
            switch Inputs.ArmatureType
                case 'external'
                    YShift = design.Rmo + design.g + (design.Ryo + 2*design.tm + 10*design.tm) + 10*Inputs.Tol;
                case 'internal'
                    YShift = design.Ryo + design.g + (design.Rbo + sum (rotorinfo.OuterWrapperThickness(2:end))) + 10*Inputs.Tol;
            end
            
            % draw the stator slots
            [FemmProblem, statorinfo] = radialfluxstator2dfemmprob ( ...
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
                'ShoeCurveControlFrac', shoecurvefrac, ...
                'NSlots', Inputs.NSlots, ...
                'YShift', YShift);
            
            [FemmProblem, statorinfo] = stator_iron_boundary (FemmProblem, design, Inputs, statorinfo, Rs, 0, YShift);
            
            [FemmProblem, statorinfo] = stator_outer_regions (FemmProblem, design, Inputs, statorinfo, Rs, GapMatInd, 0, YShift);
            
            % add all the required boundaries for the gap

            % create the boundaries
            gapboundinds = nan * ones (1, Inputs.NBoundaryPositions);
            rotorgapboundnames = cell(1, Inputs.NBoundaryPositions);
            for ind = 1:Inputs.NBoundaryPositions
                [FemmProblem, gapboundinds(ind), rotorgapboundnames{ind}] = ...
                    addboundaryprop_mfemm (FemmProblem, sprintf ('gap_bound_%d', ind), 4);
            end
            
            % shift them round by the desired amount
            statorgapboundnames = circshift (rotorgapboundnames, [0, Inputs.BoundaryShift]);
%             gapboundinds = circshift (gapboundinds, [Inputs.NBoundaryPositions, 0]);
            
            switch design.ArmatureType
                
                case 'external'
                    % obtain appropriate locations for the nodes to be added for
                    % each gap boundary
                    [n1x, n1y] = pol2cart ( linspace (0, Inputs.NPolePairs*2*design.thetap, Inputs.NBoundaryPositions+1)', ...
                                            repmat (design.Rmo + design.g/2, Inputs.NBoundaryPositions+1, 1 ) ...
                                          );
                case 'internal'
                    % obtain appropriate locations for the nodes to be added for
                    % each gap boundary
                    [n1x, n1y] = pol2cart ( linspace (0, Inputs.NPolePairs*2*design.thetap, Inputs.NBoundaryPositions+1)', ...
                                            repmat (design.Rmi - design.g/2, Inputs.NBoundaryPositions+1, 1 ) ...
                                          );
            end
            % shift the nodes to create the boundary nodes for the 
            n2x = n1x;
            n2y = n1y + YShift;

            % add the nodes for the gap boundaries
            if abs(tau-(design.thetap * 2 * Inputs.NPolePairs)) > Inputs.Tol
                [FemmProblem, ~, n1nodeids] = addnodes_mfemm ( FemmProblem, n1x, n1y);
                [FemmProblem, ~, n2nodeids] = addnodes_mfemm ( FemmProblem, n2x, n2y);
                
                % add segments with periodic boundaries on the outer parts
                [FemmProblem, ~, rotorairgapboundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Air Gap Periodic', 4);
                [FemmProblem, ~, statorairgapboundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Air Gap Periodic', 4);

                % bottom gap corner to bottom rotor corner
                [FemmProblem, commoninfo.BottomSegInds] = addsegments_mfemm (FemmProblem, ...
                                                   n1nodeids(1), ...
                                                   rotorinfo.MagnetCornerIDs(1), ...
                                                   'BoundaryMarker', rotorairgapboundname);
                
                % top gap corner to top rotor corner
                [FemmProblem, commoninfo.BottomSegInds(end+1)] = addsegments_mfemm (FemmProblem, ...
                                                   n1nodeids(end), ...
                                                   rotorinfo.MagnetCornerIDs(2), ...
                                                   'BoundaryMarker', rotorairgapboundname);
                                               
                % bottom gap corner to bottom core corner
                [FemmProblem, commoninfo.BottomSegInds] = addsegments_mfemm (FemmProblem, ...
                                                   statorinfo.node_id_set_stator_inner(1), ...
                                                   n2nodeids(1), ...
                                                   'BoundaryMarker', statorairgapboundname);
                
                % top gap corner to top core corner
                [FemmProblem, commoninfo.BottomSegInds(end+1)] = addsegments_mfemm (FemmProblem, ...
                                                   statorinfo.node_id_set_stator_inner(2), ...
                                                   n2nodeids(end), ...
                                                   'BoundaryMarker', statorairgapboundname);
                                               
                                               
            else
                [FemmProblem, ~, n1nodeids] = addnodes_mfemm ( FemmProblem, n1x(1:end-1), n1y(1:end-1));
                [FemmProblem, ~, n2nodeids] = addnodes_mfemm ( FemmProblem, n2x(1:end-1), n2y(1:end-1));
                n1nodeids(end+1) = n1nodeids(1);
                n2nodeids(end+1) = n2nodeids(1);
            end

            gaparcangle = rad2deg (Inputs.NPolePairs*2*design.thetap/Inputs.NBoundaryPositions);
            maxgapsegdegress = Inputs.NPolePairs*2*design.thetap/Inputs.NBoundaryPositions/10;
            
            % create the segments for the gap boundaries
            for ind = 1:numel (n1x)-1

               % add a new gap boundary
               
               % add the gap segments
               FemmProblem = addarcsegments_mfemm ( FemmProblem, ...
                                                    n1nodeids(ind), ...
                                                    n1nodeids(ind+1), ...
                                                    gaparcangle, ...
                                                    'MaxSegDegrees', maxgapsegdegress, ...
                                                    'BoundaryMarker', rotorgapboundnames{ind} );

               % add the gap segments
               FemmProblem = addarcsegments_mfemm ( FemmProblem, ...
                                                    n2nodeids(ind), ...
                                                    n2nodeids(ind+1), ...
                                                    gaparcangle, ...
                                                    'MaxSegDegrees', maxgapsegdegress, ...
                                                    'BoundaryMarker', statorgapboundnames{ind} );
            end
            
            
            % add the air gap labels
            switch Inputs.ArmatureType
        
                case 'external'

                    % Add block labels for the air gap
                    [labelloc(1),labelloc(2)]  = pol2cart (design.thetap, design.Rmo + design.g/4);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm ( labelloc(1,1), labelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize );
                                        
                    [labelloc(1),labelloc(2)]  = pol2cart (design.thetap, design.Rmo + 3*design.g/4);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm ( labelloc(1,1), labelloc(1,2) + YShift, ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize );

                    % add a block label for the yoke and teeth
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, Rs);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2) + YShift, ...
                                                                         yokeBlockProps);

                case 'internal'

                    % Add block labels for the air gap
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, design.Rmi-design.g/4);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize);
                                        
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, design.Rmi-3*design.g/4);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2) + YShift, ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize);

                    % add a block label for the yoke and teeth
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, Rs);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                                                         yokeBlockProps);

                otherwise
                    error('Unrecognised ArmatureType option.')

            end
            
            
        case 'Full'
            
            % we're drawing everything and theefore ignoring some inputs
            Inputs.NSlots = design.Qs;
            Inputs.NPolePairs = design.Poles / 2;
            
            % draw the stator slots
            [FemmProblem, statorinfo] = radialfluxstator2dfemmprob ( ...
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
                    'ShoeCurveControlFrac', shoecurvefrac, ...
                    'NSlots', design.Qs );
            
            [FemmProblem, statorinfo] = stator_iron_boundary (FemmProblem, design, Inputs, statorinfo, Rs, 0, 0);
            
            [FemmProblem, statorinfo] = stator_outer_regions (FemmProblem, design, Inputs, statorinfo, Rs, GapMatInd, 0, 0);
                
%             % close the last gap in the stator
%             angle = rad2deg (design.thetas-design.thetacg);
% 
%             [FemmProblem, ~] = addarcsegments_mfemm ( FemmProblem, ...
%                                                       statorinfo.OuterNodes(end), ...
%                                                       statorinfo.OuterNodes(1), ...
%                                                       rad2deg (design.thetas-design.thetacg), ...
%                                                       'MaxSegDegrees', min (angle/5, 1) );

            switch Inputs.ArmatureType
                
                case 'external'
                    
                    FemmProblem = addgroup_mfemm (FemmProblem, 'Magnet');
                    FemmProblem = addgroup_mfemm (FemmProblem, 'BackIron');
                    
%                     statorinfo.touterregion = [2*design.tm, 10*design.tm]; 
%                     
%                     % make the outer regions
%                     [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Ryo + statorinfo.touterregion(1); -(design.Ryo + statorinfo.touterregion(1))], [0; 0]);
% 
%                     [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
%                                                                        nodeids(1), ...
%                                                                        nodeids(2), ...
%                                                                        180 );
% 
%                     [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
%                                                                        nodeids(2), ...
%                                                                        nodeids(1), ...
%                                                                        180 );
%                                                                    
%                     % make the outer regions
%                     [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Ryo + statorinfo.touterregion(2); -(design.Ryo + statorinfo.touterregion(2))], [0; 0]);
% 
%                     [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
%                                                                        nodeids(1), ...
%                                                                        nodeids(2), ...
%                                                                        180 );
% 
%                     [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
%                                                                        nodeids(2), ...
%                                                                        nodeids(1), ...
%                                                                        180 );
%                                                                    
%                     % Add block labels for the outer air regions
%                     FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (design.Ryo + statorinfo.touterregion(1)/2, 0, ...
%                                             'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
%                                             'MaxArea', Inputs.OuterRegionsMeshSize(1));
% 
%                     FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (design.Ryo + statorinfo.touterregion(1) + statorinfo.touterregion(2)/2, 0, ...
%                                             'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
%                                             'MaxArea', Inputs.OuterRegionsMeshSize(2));
%                                                                    
%                     
%                     % make the outer boundary of the stator back iron
%                     [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Ryo; -design.Ryo], [0; 0]);
% 
%                     [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
%                                                                        nodeids(1), ...
%                                                                        nodeids(2), ...
%                                                                        180 );
% 
%                     [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
%                                                                        nodeids(2), ...
%                                                                        nodeids(1), ...
%                                                                        180 );

                    % make the rotor iron boundary
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Rbi; -design.Rbi], [0; 0]);

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(1), ...
                                                                       nodeids(2), ...
                                                                       180, ...
                                                                       'InGroup', FemmProblem.Groups.BackIron);

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(2), ...
                                                                       nodeids(1), ...
                                                                       180, ...
                                                                       'InGroup', FemmProblem.Groups.BackIron);
                                                                   
                    % make inner air region
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, 0.5*[design.Rbi; -design.Rbi], [0; 0]);

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(1), ...
                                                                       nodeids(2), ...
                                                                       180 );

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(2), ...
                                                                       nodeids(1), ...
                                                                       180 );
                                                                   
                    % add air label
                    FemmProblem = addblocklabel_mfemm ( FemmProblem, 0.75*design.Rbi, 0, ...
                                            'BlockType', FemmProblem.Materials(1).Name );
                                        
                    % add center label
                    FemmProblem = addblocklabel_mfemm ( FemmProblem, 0, 0, ...
                                            'BlockType', '<No Mesh>' );
                                        
                    % Add block label for the air gap
                    [labelloc(1),labelloc(2)]  = pol2cart (design.thetap, design.Rmo+design.g/2);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm ( labelloc(1,1), labelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.AirGapMeshSize );

                    % add a block label for the yoke and teeth
                    [labelloc(1),labelloc(2)] = pol2cart (design.thetap, Rs);

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                                                          yokeBlockProps);
                                                           
                case 'internal'
                    
                    
                    
                otherwise
                    
                    
            end
                
            % draw the magnets
            thistheta = 0;
            thetadiff = tau / design.Poles;
            magdirs = [0, pi];
            
            for magn = 1:design.Poles
                
                magccords = [ -design.thetam/2 + thistheta, design.Rmi;
                              design.thetam/2 + thistheta, design.Rmi;
                              -design.thetam/2 + thistheta, design.Rmo;
                              design.thetam/2 + thistheta, design.Rmo;
                              (magn-1) * design.thetap, design.Rmm ];
                          
                [magccords(:,1), magccords(:,2)] = pol2cart (magccords(:,1), magccords(:,2));
                
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, ...
                                                            magccords(1:4,1), magccords(1:4,2), ...
                                                            'InGroup', FemmProblem.Groups.Magnet);
                
                % add the magnet boundary
                [FemmProblem, inarcseg] = addarcsegments_mfemm ( FemmProblem, ...
                                                                   nodeids(1), ...
                                                                   nodeids(2), ...
                                                                   rad2deg (design.thetam), ...
                                                                   'InGroup', FemmProblem.Groups.Magnet );

                [FemmProblem, outarcseg] = addarcsegments_mfemm ( FemmProblem, ...
                                                                   nodeids(3), ...
                                                                   nodeids(4), ...
                                                                   rad2deg (design.thetam), ...
                                                                   'InGroup', FemmProblem.Groups.Magnet );
                                                               
                [FemmProblem, topsegind] = addsegments_mfemm ( FemmProblem, ...
                                                                nodeids(1), ...
                                                                nodeids(3), ...
                                                                'InGroup', FemmProblem.Groups.Magnet);
                
                [FemmProblem, botsegind] = addsegments_mfemm ( FemmProblem, ...
                                                               nodeids(2), ...
                                                               nodeids(4), ...
                                                               'InGroup', FemmProblem.Groups.Magnet);
                                                           
                if magn > 1
                    % link up the magnets
                    FemmProblem = addarcsegments_mfemm ( FemmProblem, ...
                                                         lastnodeids(2), ...
                                                         nodeids(1), ...
                                                         rad2deg (design.thetap - design.thetam), ...
                                                         'InGroup', FemmProblem.Groups.Magnet );
                    
                else
                    % add the magnet block label
                    firstnodeids = nodeids;
                end
                
                FemmProblem = addblocklabel_mfemm ( FemmProblem, ...
                                        magccords(5,1), magccords(5,2), ...
                                        'BlockType', FemmProblem.Materials(MagnetMatInd).Name, ...
                                        ... 'MaxArea', Inputs.MeshSize, ...
                                        'MagDir', rad2deg (thistheta + magdirs(1)), ...
                                        'InGroup', FemmProblem.Groups.Magnet );
                
                lastnodeids = nodeids;
                                    
                magdirs = fliplr (magdirs);
                thistheta = thistheta + thetadiff;
                
            end
            
            % final magnet link
            FemmProblem = addarcsegments_mfemm ( FemmProblem, ...
                                                 lastnodeids(2), ...
                                                 firstnodeids(1), ...
                                                 rad2deg (design.thetap - design.thetam), ...
                                                 'InGroup', FemmProblem.Groups.Magnet );
            
            % add rotor back iron label
            FemmProblem = addblocklabel_mfemm ( FemmProblem, design.Rbm, 0, ...
                                            'BlockType', FemmProblem.Materials(BackIronMatInd).Name, ...
                                            'InGroup', FemmProblem.Groups.BackIron );
                                                    
        otherwise
            
            error ('Unrecognised simulation type, valid options are ''2PoleMagnetRotation'' and ''Full''');
            
    end
    
    if Inputs.DrawCoilInsulation
        FemmProblem = addcoilinsulationlabels (FemmProblem, design, Inputs, statorinfo.InsulationLabelLocations);
    end
    
    FemmProblem = addcircuitsandcoillabels (FemmProblem, design, Inputs, CoilMatInd, statorinfo.CoilLabelLocations);
    
end

function tbboundnames = getsegbounds (FemmProblem, tbboundseginds)

    % preallocate a cell array to hold the boundary names
    tbboundnames = cell (size (tbboundseginds));
    
    for ind = 1:numel (tbboundseginds)
        
        tbboundnames{ind} = FemmProblem.Segments(tbboundseginds(ind)).BoundaryMarker;
        
    end

end


function [FemmProblem, statorinfo] = stator_iron_boundary (FemmProblem, design, Inputs, statorinfo, Rs, XShift, YShift)

    elcount = elementcount_mfemm (FemmProblem);
    
    statorirongp = getgroupnumber_mfemm (FemmProblem, 'StatorIronOutline');

    edgenodes = [];
            
    switch Inputs.ArmatureType

        case 'external'
            % single inner facing stator (magnets inside, stator outside)
            
            % create the nodes which will make up the outer region and
            % upper and lower enges of the stator
            [edgenodes(:,1), edgenodes(:,2)] = pol2cart ( ...
                                                         [ 0; ...
                                                           design.thetap * 2 * Inputs.NPolePairs; ...
                                                           0; ...
                                                           design.thetap * Inputs.NPolePairs; ...
                                                           design.thetap * 2 * Inputs.NPolePairs ...
                                                           ], ...
                                                         [ design.Rmo+design.g; ... % inner surface of stator, bottom corner
                                                           design.Rmo+design.g; ... % inner surface of stator, top corner
                                                           Rs + design.ty/2; ... % outer surface of stator, bottom node
                                                           Rs + design.ty/2; ... % outer surface of stator, middle node
                                                           Rs + design.ty/2; ... % outer surface of stator, top node
                                                           ] ...
                                                         ); 

            if abs(tau-(design.thetap * 2 * Inputs.NPolePairs)) > Inputs.Tol
                
                % add the nodes to the problem
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));

                statorinfo.node_id_set_stator_inner = [nodeids(1), nodeids(2)];
                statorinfo.node_id_set_stator_outer = [nodeids(3), nodeids(4), nodeids(5)];
                statorinfo.node_id_pair_stator_bottom_edge = [nodeids(1), nodeids(3)];
                statorinfo.node_id_pair_stator_top_edge = [nodeids(2), nodeids(5)];
            
                % add  periodic boundary for the segments on the edges of the
                % stator
                [FemmProblem, ~, statorboundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Stator Back Iron Periodic', 4);

                segprop = struct ('BoundaryMarker', statorboundname, 'InGroup', statorirongp);

                % bottom stator boundary edge segment
                [FemmProblem, statorinfo.BottomSegInds] = addsegments_mfemm (FemmProblem, ...
                                                                       statorinfo.node_id_pair_stator_bottom_edge(1), ...
                                                                       statorinfo.node_id_pair_stator_bottom_edge(2), ...
                                                                       segprop);
                % top stator boundary edge segment
                [FemmProblem, statorinfo.TopSegInds] = addsegments_mfemm (FemmProblem, ...
                                                                    statorinfo.node_id_pair_stator_top_edge(1), ...
                                                                    statorinfo.node_id_pair_stator_top_edge(2), ...
                                                                    segprop);
                                                                
                % arc from bottom slot corner to sim periodic boundary edge at bottom
                FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.node_id_set_stator_inner(1), statorinfo.OuterNodes(1), ...
                                                    rad2deg(((2*pi/design.Qs)-design.thetac(1))/2), ...
                                                    'InGroup', statorirongp);

                % arc from top slot corner to sim periodic boundary edge at bottom
                FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.OuterNodes(4), statorinfo.node_id_set_stator_inner(2), ...
                                                    rad2deg(((2*pi/design.Qs)-design.thetac(1))/2), ...
                                                    'InGroup', statorirongp);
            else
                
                % add the nodes to the problem
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes([1,2,4,5],1), edgenodes([1,2,4,5],2));

%                 statorinfo.node_id_set_stator_inner = [nodeids(1), nodeids(2), nodeids(1)];
                statorinfo.node_id_set_stator_outer = [nodeids(3), nodeids(4), nodeids(3)];
                
                % arc from first slot bottom left corner to last slot top
                % left corner
                FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.OuterNodes(4), statorinfo.OuterNodes(1), ...
                                                    rad2deg(((2*pi/design.Qs)-design.thetac(1))), ...
                                                    'InGroup', statorirongp);
                
            end

            % add arcs linking the outer segments
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                statorinfo.node_id_set_stator_outer(1), ...
                                                statorinfo.node_id_set_stator_outer(2), ...
                                                rad2deg(design.thetap*Inputs.NPolePairs));

            % put the stator iron segment in the right group
            FemmProblem.ArcSegments(end).InGroup = statorirongp;
            
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                statorinfo.node_id_set_stator_outer(2), ...
                                                statorinfo.node_id_set_stator_outer(3), ...
                                                rad2deg(design.thetap*Inputs.NPolePairs));

            % put the stator iron segment in the right group
            FemmProblem.ArcSegments(end).InGroup = statorirongp;

        case 'internal'
            % single outer facing stator (stator inside, magnets outside)

            % create the nodes which will make up the outer region and
            % upper and lower enges of the stator
            [edgenodes(:,1), edgenodes(:,2)] = pol2cart ( ...
                                                         [ 0; ...
                                                           0; ...
                                                           design.thetap * 2 * Inputs.NPolePairs; ...
                                                           design.thetap * 2 * Inputs.NPolePairs ...
                                                           ], ...
                                                         [ design.Rmi-design.g; ... % outer surface of stator, bottom corner
                                                           Rs - design.ty/2; ... % inner surface of stator, bottom node
                                                           Rs - design.ty/2; ... % inner surface of stator, top node
                                                           design.Rmi-design.g; % outer surface of stator, top corner
                                                           ] ...
                                                         ); 

            if abs(tau-(design.thetap * 2 * Inputs.NPolePairs)) > Inputs.Tol
                
                % add the nodes to the problem
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));

                statorinfo.node_id_set_stator_outer = [nodeids(1), nodeids(4)];
                statorinfo.node_id_set_stator_inner = [nodeids(2), nodeids(3)];
                statorinfo.node_id_pair_stator_bottom_edge = [nodeids(1), nodeids(2)];
                statorinfo.node_id_pair_stator_top_edge = [nodeids(4), nodeids(3)];
                
                % add  periodic boundaries for the segments on the edges of
                % the stator
                [FemmProblem, ~, statorboundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Stator Back Iron Periodic', 4);

                segprop = struct ('BoundaryMarker', statorboundname, 'InGroup', statorirongp);

                % bottom stator boundary edge segment
                [FemmProblem, statorinfo.BottomSegInds] = addsegments_mfemm (FemmProblem, ...
                                                                       statorinfo.node_id_pair_stator_bottom_edge(1), ...
                                                                       statorinfo.node_id_pair_stator_bottom_edge(2), ...
                                                                       segprop);
                % top stator boundary edge segment
                [FemmProblem, statorinfo.TopSegInds] = addsegments_mfemm (FemmProblem, ...
                                                                    statorinfo.node_id_pair_stator_top_edge(1), ...
                                                                    statorinfo.node_id_pair_stator_top_edge(2), ...
                                                                    segprop);

                % arc from bottom slot corner to sim periodic boundary edge at bottom
                FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.node_id_set_stator_inner(1), statorinfo.OuterNodes(1), ...
                                                    rad2deg(((2*pi/design.Qs)-design.thetac(1))/2), ...
                                                    'InGroup', statorirongp);

                % arc from top slot corner to sim periodic boundary edge at bottom
                FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.OuterNodes(4), statorinfo.node_id_set_stator_inner(2), ...
                                                    rad2deg(((2*pi/design.Qs)-design.thetac(1))/2), ...
                                                    'InGroup', statorirongp);
                                            
            else
                
                % add the nodes to the problem
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));

                statorinfo.node_id_set_stator_outer = [nodeids(1), nodeids(4)];
                statorinfo.node_id_set_stator_inner = [nodeids(2), nodeids(3)];
                statorinfo.node_id_pair_stator_bottom_edge = [nodeids(1), nodeids(2)];
                statorinfo.node_id_pair_stator_top_edge = [nodeids(4), nodeids(3)];
                
                % arc from first slot bottom left corner to last slot top
                % left corner
                FemmProblem = addarcsegments_mfemm (FemmProblem, statorinfo.OuterNodes(4), statorinfo.OuterNodes(1), ...
                                                    rad2deg(((2*pi/design.Qs)-design.thetac(1))), ...
                                                    'InGroup', statorirongp);
                
            end

            % add arcs linking the outer segments
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                statorinfo.node_id_set_stator_outer(1), ...
                                                statorinfo.node_id_set_stator_outer(2), ...
                                                rad2deg(design.thetap*Inputs.NPolePairs));

            % put the stator iron segment in the right group
            FemmProblem.ArcSegments(end).InGroup = statorirongp;
            
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                statorinfo.node_id_set_stator_outer(1), ...
                                                statorinfo.node_id_set_stator_outer(2), ...
                                                rad2deg(design.thetap*Inputs.NPolePairs));

            % put the stator iron segment in the right group
            FemmProblem.ArcSegments(end).InGroup = statorirongp;



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
    
    FemmProblem = translatenewelements_mfemm (FemmProblem, elcount, XShift, YShift);

end


function [FemmProblem, statorinfo] = stator_outer_regions (FemmProblem, design, Inputs, statorinfo, Rs, GapMatInd, XShift, YShift)

    elcount = elementcount_mfemm (FemmProblem);
    
    edgenodes = [];
            
    switch Inputs.ArmatureType

        case 'external'
            % single inner facing stator (magnets inside, stator outside)

            % the sizes of the stator outer air regions
            if isempty (Inputs.StatorOuterRegionSize)
                Inputs.StatorOuterRegionSize = [2*design.tm, 10*design.tm];                
            end
            if isempty (Inputs.StatorOuterRegionsMeshSize)
                Inputs.StatorOuterRegionsMeshSize = [ choosemesharea_mfemm(design.tm, (Rs*design.thetap), 1/5), ...
                                                      repmat(-1,1,numel(Inputs.StatorOuterRegionSize)-1) ];
            end
            assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionsMeshSize), ...
                'RENEWNET:slottedfemmproble_radial:nstatormeshsizes', ...
                'Number of supplied stator outer region mesh sizes does not match number of outer region sizes.');
            
            if isempty (Inputs.StatorOuterRegionMaterials)
                % FemmProblem.Materials(GapMatInd).Name
                Inputs.StatorOuterRegionMaterials = repmat ({FemmProblem.Materials(GapMatInd).Name}, 1, numel(Inputs.StatorOuterRegionSize));
            end
            assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionMaterials), ...
                'RENEWNET:slottedfemmproble_radial:nstatormaterials', ...
                'Number of supplied stator outer region material names does not match number of outer region sizes.');
                    
            statorinfo.touterregion = Inputs.StatorOuterRegionSize;
            
            % calculate the radial positions of the outer region boundaries
            statorinfo.routerregion = cumsum ([Rs + design.ty/2, statorinfo.touterregion]);
            
            % get a suitable tolerance to determine if we are making a full
            % circle or not
            FullCircleAngleTol = Inputs.Tol / statorinfo.routerregion(1);
            
            sectionangle = design.thetap * 2 * Inputs.NPolePairs;
            
            statorinfo.node_id_sets_stator_outer_region = statorinfo.node_id_set_stator_outer;
            
            for ind = 1:numel (statorinfo.touterregion)

                if (sectionangle < tau()) && ((tau()-sectionangle) > FullCircleAngleTol)
                    % there is not a full circle, so we will be making
                    % boundaries at the top and bottom. Test above is
                    % checking angle and tolerance based on actual
                    % distance, not angle
                    
                    [edgenodes(:,1), edgenodes(:,2)] = ...
                        pol2cart ( [ 0; design.thetap * Inputs.NPolePairs; design.thetap * 2 * Inputs.NPolePairs; ], ...
                                     repmat (statorinfo.routerregion(ind+1), [3,1] ) );

                    % add the appropriate nodes to the problem
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));

                    statorinfo.node_id_sets_stator_outer_region = [statorinfo.node_id_sets_stator_outer_region;
                                                                   nodeids];
                                                               
%                     statorinfo.node_id_set_stator_outer_region_2 = [nodeids(4), nodeids(5), nodeids(6)];

                    % add segments with periodic boundaries on the outer parts
                    [FemmProblem, ~, boundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Stator Outer Periodic', 4);
%                     [FemmProblem, ~, region2boundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Stator Outer Periodic', 4);

                    segprops = struct ('BoundaryMarker', boundname, 'InGroup', 0);

                    % bottom seg
                    [FemmProblem, statorinfo.BottomSegInds(ind)] ...
                        = addsegments_mfemm ( FemmProblem, ...
                                              statorinfo.node_id_sets_stator_outer_region(ind,1), ...
                                              statorinfo.node_id_sets_stator_outer_region(ind+1,1), ...
                                              segprops );
                    % top seg
                    [FemmProblem, statorinfo.TopSegInds(ind)] ...
                        = addsegments_mfemm ( FemmProblem, ...
                                              statorinfo.node_id_sets_stator_outer_region(ind,3), ...
                                              statorinfo.node_id_sets_stator_outer_region(ind+1,3), ...
                                              segprops );

                elseif (sectionangle > tau()) && ((sectionangle - tau()) > FullCircleAngleTol)
                    error ('RENEWNET:slottedfemmprob_radial:toomanypoles', ...
                        'You are attempting to draw too many poles, angle is greater than a complete circle.')
                else
                    % there is a full circle, so we'll be linking up the top
                    % and bottom to itelf
                    [edgenodes(:,1), edgenodes(:,2)] = ...
                        pol2cart ( [ 0; design.thetap * Inputs.NPolePairs; ], ...
                                     repmat (statorinfo.routerregion(ind+1), [2,1] ) );
                    
                    % add the nodes to the problem
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));

                    statorinfo.node_id_sets_stator_outer_region ...
                        = [statorinfo.node_id_sets_stator_outer_region;
                           nodeids(1), nodeids(2), nodeids(1)];
                       
%                     statorinfo.node_id_set_stator_outer_region_2 = [nodeids(3), nodeids(4), nodeids(3)];

                end

                % add arcs linking the outer segments
                FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                   [ statorinfo.node_id_sets_stator_outer_region(ind+1,1), ...
                                                     statorinfo.node_id_sets_stator_outer_region(ind+1,2), ...
                                                     ], ...
                                                   [ statorinfo.node_id_sets_stator_outer_region(ind+1,2), ...
                                                     statorinfo.node_id_sets_stator_outer_region(ind+1,3), ...
                                                     ], ...
                                                   rad2deg(repmat(design.thetap*Inputs.NPolePairs,1,2)));

                % Add block labels for the outer air region
                [labelloc(1),labelloc(2)] = pol2cart (design.thetap * Inputs.NPolePairs, mean(statorinfo.routerregion(ind:ind+1)));

                
                FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', Inputs.StatorOuterRegionMaterials{ind}, ...
                                        'MaxArea', Inputs.StatorOuterRegionsMeshSize(ind));
% 
%                 [labelloc(1),labelloc(2)]  = pol2cart(design.thetap * Inputs.NPolePairs, Rs + design.ty/2 + statorinfo.touterregion(1) + statorinfo.touterregion(2)/2);
% 
%                 FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
%                                         'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
%                                         'MaxArea', Inputs.OuterRegionsMeshSize(2));

            end
            
        case 'internal'
            % single outer facing stator (stator inside, magnets outside)

            statorinfo.touterregion = [0.5*design.Ryi, 0.8*design.Ryi]; 
                                                     
            [edgenodes(:,1), edgenodes(:,2)] = pol2cart ( ...
                                                         [ 0; ...
                                                           design.thetap * Inputs.NPolePairs; ...
                                                           design.thetap * 2 * Inputs.NPolePairs; ...
                                                           0; ...
                                                           design.thetap * Inputs.NPolePairs; ...
                                                           design.thetap * 2 * Inputs.NPolePairs; ...
                                                           ], ...
                                                         [ statorinfo.touterregion(1); ... % outer surface of first inner region, bottom node
                                                           statorinfo.touterregion(1); ... % outer surface of first inner region, middle node
                                                           statorinfo.touterregion(1); ... % outer surface of first inner region, top node
                                                           statorinfo.touterregion(2); ... % outer surface of second outer region, bottom node
                                                           statorinfo.touterregion(2); ... % outer surface of second outer region, middle node
                                                           statorinfo.touterregion(2); ... % outer surface of second outer region, top node
                                                           ] ...
                                                         ); 



            if abs(tau-(design.thetap * 2 * Inputs.NPolePairs)) > Inputs.Tol
                
                % add the appropriate nodes to the problem
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes(:,1), edgenodes(:,2));

                statorinfo.node_id_set_stator_outer_region_1 = [nodeids(1), nodeids(2), nodeids(3)];
                statorinfo.node_id_set_stator_outer_region_2 = [nodeids(4), nodeids(5), nodeids(6)];

                % add segments with periodic boundaries on the outer parts
                [FemmProblem, ~, region1boundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Stator Outer Periodic', 4);
                [FemmProblem, ~, region2boundname] = addboundaryprop_mfemm (FemmProblem, 'Radial Stator Outer Periodic', 4);

                segprops = struct ('BoundaryMarker', {region1boundname, region2boundname}, 'InGroup', {0, 0});
                
                % bottom segs
                [FemmProblem, statorinfo.BottomSegInds(end+1:end+2)] = addsegments_mfemm (FemmProblem, ...
                                                [statorinfo.node_id_set_stator_outer(1), statorinfo.node_id_set_stator_outer_region_1(1)], ...
                                                [statorinfo.node_id_set_stator_outer_region_1(1) statorinfo.node_id_set_stator_outer_region_2(1)], ...
                                                segprops);
                % top segs
                [FemmProblem, statorinfo.TopSegInds(end+1:end+2)] = addsegments_mfemm (FemmProblem, ...
                                                [statorinfo.node_id_set_stator_outer(3), statorinfo.node_id_set_stator_outer_region_1(3)], ...
                                                [statorinfo.node_id_set_stator_outer_region_1(3) statorinfo.node_id_set_stator_outer_region_2(3)], ...
                                                segprops);

            else
                
                % add the nodes to the problem
                [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, edgenodes([1,2,4,5],1), edgenodes([1,2,4,5],2));

                statorinfo.node_id_set_stator_outer_region_1 = [nodeids(1), nodeids(2), nodeids(1)];
                statorinfo.node_id_set_stator_outer_region_2 = [nodeids(3), nodeids(4), nodeids(3)];
                
            end

            % add arcs linking the outer segments
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                               [ statorinfo.node_id_set_stator_outer_region_1(1), ...
                                                 statorinfo.node_id_set_stator_outer_region_1(2), ...
                                                 statorinfo.node_id_set_stator_outer_region_2(1), ...
                                                 statorinfo.node_id_set_stator_outer_region_2(2), ...
                                                 ], ...
                                               [ statorinfo.node_id_set_stator_outer_region_1(2), ...
                                                 statorinfo.node_id_set_stator_outer_region_1(3), ...
                                                 statorinfo.node_id_set_stator_outer_region_2(2), ...
                                                 statorinfo.node_id_set_stator_outer_region_2(3)], ...
                                               rad2deg(repmat(design.thetap*Inputs.NPolePairs,1,4)));


            % Add block labels for the outer air regions
            [labelloc(1),labelloc(2)]  = pol2cart (design.thetap * Inputs.NPolePairs, ...
                                            statorinfo.touterregion(1) + (Rs - design.ty/2 - statorinfo.touterregion(1))/2 );

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                    'MaxArea', Inputs.OuterRegionsMeshSize(1));

            [labelloc(1),labelloc(2)]  = pol2cart (design.thetap * Inputs.NPolePairs, ...
                                                    statorinfo.touterregion(2) + (statorinfo.touterregion(1) - statorinfo.touterregion(2))/2);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                    'MaxArea', Inputs.OuterRegionsMeshSize(2));

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
    
    % move the new nodes and labels
    FemmProblem = translatenewelements_mfemm (FemmProblem, elcount, XShift, YShift);
            
end


% function FemmProblem = shift_new_nodes_and_labels (FemmProblem, oldelcount, XShift, YShift)
% % shift newly added nodes blocks and labels
% %
% % Syntax
% %
% % FemmProblem = shiftnewelements_mfemm (FemmProblem, oldelcount, XShift, YShift)
% %
% % Description
% %
% % shiftnewelements_mfemm moves all nodes and block lables added since the
% % element count provided in oldelcount was performed by the specified
% % amount in the x and y direction
% %
% % 
% 
%     if XShift ~= 0 || YShift ~= 0
%         
%         newelcount = elementcount_mfemm (FemmProblem);
%         
%         nodeids = (oldelcount.NNodes):(newelcount.NNodes-1);
%         
%         FemmProblem = translatenodes_mfemm(FemmProblem, XShift, YShift, nodeids);
%         
%         blockids = (oldelcount.NBlockLabels):(newelcount.NBlockLabels-1);
%                  
%         FemmProblem = translateblocklabels_mfemm(FemmProblem, XShift, YShift, blockids);
%         
%     end
% 
% end


function FemmProblem = addcoilinsulationlabels (FemmProblem, design, Inputs, inslabellocs)
    
   if strncmpi (Inputs.SimType, 'm', 1)
       [FemmProblem, matinds] = addmaterials_mfemm (FemmProblem, ...
                                                    {design.MagFEASimMaterials.CoilInsulation}, ...
                                                    'MaterialsLibrary', Inputs.MaterialsLibrary);
   elseif strncmpi (Inputs.SimType, 'h', 1)
       [FemmProblem, matinds] = addmaterials_mfemm (FemmProblem, ...
                                                    {design.HeatFEASimMaterials.CoilInsulation}, ...
                                                    'MaterialsLibrary', Inputs.MaterialsLibrary );
   end

    % add the coil insulation block labels
    coilinsBlockProps.BlockType = FemmProblem.Materials(matinds).Name;
    coilinsBlockProps.MaxArea = Inputs.CoilInsRegionMeshSize;
    coilinsBlockProps.InCircuit = '';
    coilinsBlockProps.InGroup = Inputs.CoilGroup;

    % add insulation labels if requested

    for indi = 1:size (inslabellocs)
        FemmProblem = addblocklabel_mfemm( FemmProblem, ...
                                           inslabellocs(indi,1), ...
                                           inslabellocs(indi,2), ...
                                           coilinsBlockProps );
    end
        

    
end


function FemmProblem = addcircuitsandcoillabels (FemmProblem, design, Inputs, coilmatind, coillabellocs)

    % add circuits for each winding phase
    for i = 1:design.Phases
        cname = num2str(i);
        if ~hascircuit_mfemm (FemmProblem, cname)
            FemmProblem = addcircuit_mfemm (FemmProblem, cname);
        end
        FemmProblem = setcircuitcurrent (FemmProblem, cname, Inputs.CoilCurrent(i));
    end
    
    coilBlockProps.BlockType = FemmProblem.Materials(coilmatind).Name;
    coilBlockProps.MaxArea = Inputs.CoilRegionMeshSize;
    coilBlockProps.InCircuit = '';
    coilBlockProps.InGroup = Inputs.CoilGroup;

    % draw the positive part of the coil circuit
    coilBlockProps.Turns = design.CoilTurns;
    
    % add block labels for the coils
    row = 1;
    for slotn = 1:Inputs.NSlots

        for layern = 1:Inputs.NWindingLayers

            coilBlockProps.InCircuit = num2str(abs(design.WindingLayout.Phases(slotn,layern)));
            coilBlockProps.Turns = design.CoilTurns * sign (design.WindingLayout.Phases(slotn,layern));
            
            FemmProblem = addblocklabel_mfemm( FemmProblem, ...
                                               coillabellocs(row,1), ...
                                               coillabellocs(row,2), ...
                                               coilBlockProps);

            row = row + 1;

        end

    end
    
end
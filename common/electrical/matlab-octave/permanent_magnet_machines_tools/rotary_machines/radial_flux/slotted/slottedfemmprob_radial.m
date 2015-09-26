function [FemmProblem, coillabellocs] = slottedfemmprob_radial(design, varargin)
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
%   'BoundaryPositions' = 1;
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
    Inputs.BoundaryPositions = 1;
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
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm (design.g, (design.Rmm*design.thetap), 1/10);
    Inputs.DrawOuterRegions = true;
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
    Inputs.NSlots = Inputs.NPolePairs*2*design.Qs/design.Poles;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
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
    coillabellocs = [];
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
    
    switch Inputs.DrawingType
        
        case 'MagnetRotation'
    
            for ind = 1:Inputs.NPolePairs

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
                    'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
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
                    'ArmatureIronGroup', Inputs.ArmatureBackIronGroup, ...
                    'Tol', Inputs.Tol, ...
                    'DrawCoilInsulation', Inputs.DrawCoilInsulation, ...
                    'CoilInsulationThickness', design.CoilInsulationThickness, ...
                    'CoilBaseFraction', coilbasefrac, ...
                    'ShoeCurveControlFrac', shoecurvefrac, ...
                    'NSlots', Inputs.NSlots, ...
                    'YShift', YShift );

                coillabellocs = [coillabellocs; statorinfo.CoilLabelLocations];

                % Complete the design using the common radial drawing function
                Inputs.AddAllCoilsToCircuits = true;
                Inputs.StartSlot = lastslot;
                [FemmProblem, commoninfo] = slottedcommonfemmprob_radial ( FemmProblem, ...
                                                            design, ...
                                                            Inputs, ...
                                                            rotorinfo.MagnetCornerIDs, ... magcornerids, ...
                                                            Rs, ...
                                                            statorinfo.CoilLabelLocations, ...
                                                            statorinfo.InsulationLabelLocations, ...
                                                            statorinfo.OuterNodes, ...
                                                            design.thetap, ...
                                                            BackIronMatInd, ...
                                                            YokeMatInd, ...
                                                            CoilMatInd, ...
                                                            GapMatInd, ...
                                                            rotorinfo.LinkTopBottom, ...
                                                            XShift, ...
                                                            YShift );



                lastslot = lastslot + Inputs.NSlots;
                YShift = YShift + (2.01 * outerR);
                XShift = 0;

                tbboundseginds = [ [rotorinfo.TopSegInds, commoninfo.TopSegInds];
                                   [rotorinfo.BottomSegInds, commoninfo.BottomSegInds];
                                   tbboundseginds ];

            end

            % rearrange the boundaries to link everything up correctly

            % first get all the boundary IDs in the appropriate segments
            tbboundnames = getsegbounds (FemmProblem, tbboundseginds);

            % shift them round by one
            tbboundnames = circshift (tbboundnames, [-1, 0]);

            % replace the boundaries
            for ind = 1:numel (tbboundnames)

                FemmProblem.Segments(tbboundseginds(ind)).BoundaryMarker = tbboundnames{ind};

            end
            
            
        case 'Full'
            
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
                
            % close the last gap in the stator
            angle = rad2deg (design.thetas-design.thetacg);

            [FemmProblem, ~] = addarcsegments_mfemm ( FemmProblem, ...
                                                      statorinfo.OuterNodes(end), ...
                                                      statorinfo.OuterNodes(1), ...
                                                      rad2deg (design.thetas-design.thetacg), ...
                                                      'MaxSegDegrees', min (angle/5, 1) );

            switch Inputs.ArmatureType
                
                case 'external'
                    
                    FemmProblem = addgroup_mfemm (FemmProblem, 'Magnet');
                    FemmProblem = addgroup_mfemm (FemmProblem, 'BackIron');
                    
                    routerregion = [2*design.tm, 10*design.tm]; 
                    
                    % make the outer regions
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Ryo + routerregion(1); -(design.Ryo + routerregion(1))], [0; 0]);

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(1), ...
                                                                       nodeids(2), ...
                                                                       180 );

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(2), ...
                                                                       nodeids(1), ...
                                                                       180 );
                                                                   
                    % make the outer regions
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Ryo + routerregion(2); -(design.Ryo + routerregion(2))], [0; 0]);

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(1), ...
                                                                       nodeids(2), ...
                                                                       180 );

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(2), ...
                                                                       nodeids(1), ...
                                                                       180 );
                                                                   
                    % Add block labels for the outer air regions
                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (design.Ryo + routerregion(1)/2, 0, ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.OuterRegionsMeshSize(1));

                    FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm (design.Ryo + routerregion(1) + routerregion(2)/2, 0, ...
                                            'BlockType', FemmProblem.Materials(GapMatInd).Name, ...
                                            'MaxArea', Inputs.OuterRegionsMeshSize(2));
                                                                   
                    
                    % make the outer boundary of the rotor back iron
                    [FemmProblem, ~, nodeids] = addnodes_mfemm (FemmProblem, [design.Ryo; -design.Ryo], [0; 0]);

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(1), ...
                                                                       nodeids(2), ...
                                                                       180 );

                    [FemmProblem, arcseginds] = addarcsegments_mfemm ( FemmProblem, ...
                                                                       nodeids(2), ...
                                                                       nodeids(1), ...
                                                                       180 );

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
                                                           
            % Complete the design using the common radial drawing function
            Inputs.AddAllCoilsToCircuits = true;
            Inputs.StartSlot = lastslot;
            Inputs.NSlots = design.Qs;
            FemmProblem = slottedcommonfemmprob_radial ( FemmProblem, ...
                                                        design, ...
                                                        Inputs, ...
                                                        [], ...
                                                        Rs, ...
                                                        statorinfo.CoilLabelLocations, ...
                                                        statorinfo.InsulationLabelLocations, ...
                                                        statorinfo.OuterNodes, ...
                                                        design.thetap, ...
                                                        BackIronMatInd, ...
                                                        YokeMatInd, ...
                                                        CoilMatInd, ...
                                                        GapMatInd, ...
                                                        [], ...
                                                        0, ...
                                                        0 );
                                                    
        case 'LinkedGapBoundary'
                                                        

        otherwise
            
            error ('Unrecognised simulation type, valid options are ''2PoleMagnetRotation'' and ''Full''');
            
    end
    
    
end

function tbboundnames = getsegbounds (FemmProblem, tbboundseginds)

    % preallocate a cell array to hold the boundary names
    tbboundnames = cell (size (tbboundseginds));
    
    for ind = 1:numel (tbboundseginds)
        
        tbboundnames{ind} = FemmProblem.Segments(tbboundseginds(ind)).BoundaryMarker;
        
    end

end

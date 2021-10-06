function [FemmProblem, outermagsep, coillabellocs, yokenodeids] = slottedfemmprob_torus(design, varargin)
% creates a FemmProblem structure for a slotted torus axial flux permanent
% magnet machine

    magsep = 2*design.g + 2*design.tc + 2*design.tsb + design.ty;
    
    Inputs.NStages = 1;
    Inputs.NWindingLayers = 1;
    Inputs.CoilCurrent = 0;
    Inputs.MagArrangement = 'NN';
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.Rmo - design.Rmi);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup = [];
    Inputs.BackIronGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm(design.tm, design.taumm, 1/40);
    Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(min(design.tbi), design.taupm, 1/40);
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, design.taumm, 1/10), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm(magsep, design.taupm, 1/50);
    Inputs.ShoeGapRegionMeshSize = choosemesharea_mfemm(design.tsg, design.tausgm, 1/50);
    Inputs.YokeRegionMeshSize = min( choosemesharea_mfemm(design.ty, 2*design.taupm, 1/40), ...
                                      choosemesharea_mfemm(design.tc, design.taucs, 1/40)  );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm(design.tc, design.taucs);
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(design.taupm, Inputs.Position, Inputs.FractionalPolePosition, Inputs.RotorAnglePosition);
    
    Inputs.NSlots = 2*design.Qs/design.Poles;
    
    FemmProblem = Inputs.FemmProblem;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    [FemmProblem, matinds] = addmaterials_mfemm(FemmProblem, ...
        {design.MagFEASimMaterials.Magnet, design.MagFEASimMaterials.FieldBackIron, design.MagFEASimMaterials.ArmatureYoke, design.MagFEASimMaterials.ArmatureCoil});
                 
    MagnetMatInd = matinds(1);
    BackIronMatInd = matinds(2);
    YokeMatInd = matinds(3);
    CoilMatInd = matinds(4);
    
    % draw the torus rotor according to the spec in the design strucure
    [FemmProblem, outermagsep, innerstagewidth] = torusrotor2dfemmprob( ...
        design.taupm, design.taumm, design.tm, design.tbi, magsep, ...
        'NStages', Inputs.NStages, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', MagnetMatInd, ...
        'BackIronMaterial', BackIronMatInd, ...
        'OuterRegionsMaterial', 1, ... % Air
        'MagnetSpaceMaterial', 1, ... % Air
        'MagnetGroup', Inputs.MagnetGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'BackIronGroup', Inputs.BackIronGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'Tol', Inputs.Tol);
    
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
        'Tol', Inputs.Tol);
    
    % link the rotor stages along the top and bottom, add antiperiodic
    % boundaries, and add the coil regions.
    gapedgenodes = [-outermagsep/2, 0;
                    -outermagsep/2+magsep, 0;
                    -outermagsep/2+magsep, 2*design.taupm;
                    -outermagsep/2, 2*design.taupm];
                
    FemmProblem = slottedcommonfemmprob_torus(FemmProblem, design, ...
        Inputs, magsep, gapedgenodes, innerstagewidth, coillabellocs, ...
        yokenodeids, design.taupm, BackIronMatInd, YokeMatInd, CoilMatInd);
    
%                 
%     % define the block properties of the core region
%     yokeBlockProps.BlockType = FemmProblem.Materials(BackIronMatInd).Name;
%     yokeBlockProps.MaxArea = Inputs.BackIronRegionMeshSize;
%     yokeBlockProps.InCircuit = '';
%     yokeBlockProps.InGroup = 0;
%         
%     % Prototype an array of segprops structures
%     SegProps.MaxSideLength = -1;
%     SegProps.Hidden = 0;
%     SegProps.InGroup = 0;
%     SegProps.BoundaryMarker = '';
%     
%     SegProps = repmat(SegProps, 1, 4);
%     
%     coilBlockProps.BlockType = FemmProblem.Materials(CoilMatInd).Name;
%     coilBlockProps.MaxArea = Inputs.CoilRegionMeshSize;
%     coilBlockProps.InCircuit = '';
%     coilBlockProps.InGroup = Inputs.CoilGroup;
% 
%     % draw the positive part of the coil circuit
%     coilBlockProps.Turns = design.CoilTurns;
%     
% %     corexpos = -outermagsep/2 + design.g + design.tc;
% 
%     for i = 1:Inputs.NStages
%         
%         % get the node ids of the air gap corners on the rotor stages
%         gapcornernodeids = findnode_mfemm(FemmProblem, gapedgenodes);
%         
%         % add four new nodes at the boundary of the air gap and the teeth
%         % on the armature
%         newgapnodes = [gapedgenodes(1,1) + design.g, gapedgenodes(1,2);
%                        gapedgenodes(2,1) - design.g, gapedgenodes(2,2);
%                        gapedgenodes(3,1) - design.g, gapedgenodes(3,2)
%                        gapedgenodes(4,1) + design.g, gapedgenodes(4,2)];
%                    
%         [FemmProblem, outeryokenodeinds, outeryokenodeids] = ...
%             addnodes_mfemm(FemmProblem, newgapnodes(:,1), newgapnodes(:,2));
%         
%         % add a three new periodic boundaries for the top and bottom of the
%         % air gap and core regions
%         [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Stator Air Gap Periodic', 4);
%         [FemmProblem, boundind(2)] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Stator Back Iron Periodic', 4);
%         [FemmProblem, boundind(3)] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Stator Air Gap Periodic', 4);
%         
%         % add the segments for the top and bottom of the core       
%         [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, ...
%                                     [outeryokenodeids(1); outeryokenodeids(4)], ...
%                                     [outeryokenodeids(2); outeryokenodeids(3)], ...
%                                     'BoundaryMarker', FemmProblem.BoundaryProps(boundind(2)).Name);
%                 
%         % join up the air gaps at the top and bottom
%         
%         % bottom left gap corner to bottom left core corner
%         FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(1), outeryokenodeids(1), ...
%                                 'BoundaryMarker', FemmProblem.BoundaryProps(boundind(1)).Name);
% 
%         % top left gap corner to top left core corner
%         FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(4), outeryokenodeids(4), ...
%                                 'BoundaryMarker', FemmProblem.BoundaryProps(boundind(1)).Name);   
% 
%         % bottom right gap corner to bottom right core corner
%         FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(2), outeryokenodeids(2), ...
%                                 'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);
% 
%         % top right gap corner to top right core corner
%         FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(3), outeryokenodeids(3), ...
%                                 'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);  
% 
%         % close the tooth and air boundaries to complete the core
%         FemmProblem = addsegments_mfemm(FemmProblem, yokenodeids(i,1:4), outeryokenodeids(1:4));  
%                             
% %         corexpos = corexpos + innerstagewidth;
%         
%         % Add block labels for the air gaps
%         labelloc = [gapedgenodes(1,1) + design.g/2, design.taupm];
% 
%         FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
%                                         'BlockType', FemmProblem.Materials(1).Name, ...
%                                         'MaxArea', Inputs.AirGapMeshSize);
%                                     
%         labelloc = [gapedgenodes(2,1) - design.g/2, design.taupm];
% 
%         FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
%                                         'BlockType', FemmProblem.Materials(1).Name, ...
%                                         'MaxArea', Inputs.AirGapMeshSize);
%                                     
%         % add a block label for the yoke and teeth
%         labelloc = [gapedgenodes(1,1) + magsep/2, design.taupm];
% 
%         FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
%                                         'BlockType', FemmProblem.Materials(YokeMatInd).Name, ...
%                                         'MaxArea', Inputs.YokeRegionMeshSize);
%         
%         % add block labels for the coils
%         j = 2*(i-1) + 1;
% 
%         row = 1;
% 
%         for k = 1:2*slotsperpole
% 
%             for n = 1:Inputs.NWindingLayers
% 
%                 FemmProblem = addblocklabel_mfemm(FemmProblem, ...
%                     coillabellocs(row,j), coillabellocs(row,j+1), ...
%                     coilBlockProps);
% 
%                 row = row + 1;
% 
%                 FemmProblem = addblocklabel_mfemm(FemmProblem, ...
%                     coillabellocs(row,j), coillabellocs(row,j+1), ...
%                     coilBlockProps);
% 
%                 row = row + 1;
% 
%             end
% 
%         end
%         
%         % shift the nodes to the next location
%         gapedgenodes(:,1) = gapedgenodes(:,1) + innerstagewidth;
%         
%     end
%     
%     % add circuits for each phase
%     for i = 1:design.Phases
%        FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i));
%     end

end

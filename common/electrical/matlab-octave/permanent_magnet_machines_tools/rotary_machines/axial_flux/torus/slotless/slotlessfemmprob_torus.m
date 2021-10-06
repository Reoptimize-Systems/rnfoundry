function [FemmProblem, outermagsep] = slotlessfemmprob_torus(design, varargin)
% creates a FemmProblem structure for a slotless torus axial flux permanent
% magnet machine

    magsep = 2*design.g + 2*design.tc + design.ty;
    
    Inputs.NStages = 1;
    Inputs.DrawCoils = true;
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
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm(design.tm, design.taumm, 1/30);
    Inputs.BackIronRegionMeshSize = [];
    Inputs.BackIronRegionMeshFactor = 30;
    Inputs.OuterRegionsMeshSize = [choosemesharea_mfemm(design.tm, design.taumm, 1/10), -1];
    Inputs.AirGapMeshSize = choosemesharea_mfemm(magsep, design.taupm, 1/30);
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.BackIronRegionMeshSize)
        if numel(design.tbi) > 1
            if any(design.tbi < Inputs.Tol)
                Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(max(design.tbi), design.taupm, 1/Inputs.BackIronRegionMeshFactor);
            else
                Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(min(design.tbi), design.taupm, 1/Inputs.BackIronRegionMeshFactor);
            end
        else
            Inputs.BackIronRegionMeshSize = choosemesharea_mfemm(design.tbi, design.taupm, 1/Inputs.BackIronRegionMeshFactor);
        end
    end
    
    % Get the planar position from the position specification
    Inputs.Position = planarrotorpos(design.taupm, Inputs.Position, Inputs.FractionalPolePosition, Inputs.RotorAnglePosition);
    
    FemmProblem = Inputs.FemmProblem;
    
    elcount = elementcount_mfemm(FemmProblem);
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    FemmProblem.Materials = [FemmProblem.Materials, ...
                             matstr2matstruct_mfemm( {design.MagFEASimMaterials.Magnet, ...
                                                      design.MagFEASimMaterials.FieldBackIron} )];
    
    BackironMatInd = elcount.NMaterials + 2;
    
    % draw the torus rotor according to the spec in the design strucure
    [FemmProblem, outermagsep, innerstagewidth] = torusrotor2dfemmprob( ...
        design.taupm, design.taumm, design.tm, design.tbi, magsep, ...
        'NStages', Inputs.NStages, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', numel(FemmProblem.Materials)-1, ...
        'BackIronMaterial', numel(FemmProblem.Materials), ...
        'OuterRegionsMaterial', 1, ... % Air
        'MagnetSpaceMaterial', 1, ... % Air
        'MagnetGroup', Inputs.MagnetGroup, ...
        'BackIronGroup', Inputs.BackIronGroup, ...
        'MagnetSpaceGroup', Inputs.MagnetSpaceGroup, ...
        'MagnetRegionMeshSize', Inputs.MagnetRegionMeshSize, ...
        'BackIronRegionMeshSize', Inputs.BackIronRegionMeshSize, ...
        'OuterRegionsMeshSize', Inputs.OuterRegionsMeshSize, ...
        'Position', Inputs.Position, ...
        'Tol', Inputs.Tol);
    
    % link the rotor stages along the top and bottom, and add antiperiodic
    % boundaries

    
    gapedgenodes = [-outermagsep/2, 0;
                    -outermagsep/2+magsep, 0;
                    -outermagsep/2, 2*design.taupm;
                    -outermagsep/2+magsep, 2*design.taupm];
                
    % define the block properties of the core region
    yokeBlockProps.BlockType = FemmProblem.Materials(BackironMatInd).Name;
    yokeBlockProps.MaxArea = Inputs.BackIronRegionMeshSize;
    yokeBlockProps.InCircuit = '';
    yokeBlockProps.InGroup = 0;
        
    % Prototype an array of segprops structures
    SegProps.MaxSideLength = -1;
    SegProps.Hidden = 0;
    SegProps.InGroup = 0;
    SegProps.BoundaryMarker = '';
    
    SegProps = repmat(SegProps, 1, 4);
    
    corexpos = -outermagsep/2 + design.g + design.tc;

    for i = 1:Inputs.NStages
        
        % get the node ids of the air gap corners
        gapcornernodeids = findnode_mfemm(FemmProblem, gapedgenodes);
        
        % add a three new periodic boundaries for the top and bottom of the
        % air gap and core regions
        [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Rect Mags Air Gap Periodic', 4);
        [FemmProblem, boundind(2)] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Rect Mags Air Gap Periodic', 4);
        [FemmProblem, boundind(3)] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Rect Mags Air Gap Periodic', 4);

        % Set up the segment properties for the top and bottom of the core
        SegProps(1).BoundaryMarker = FemmProblem.BoundaryProps(boundind(2)).Name;
        SegProps(3).BoundaryMarker = SegProps(1).BoundaryMarker;
        
        % Now draw the coils if requested
        if i == 1 && Inputs.DrawCoils
            
            % create node coordinates for the yoke going clockwise around
            % the yoke corners starting from the bottom left 
            yokecoords = [-outermagsep/2 + design.g + design.tc, 0;
                          -outermagsep/2 + design.g + design.tc + design.ty, 0;
                          -outermagsep/2 + design.g + design.tc + design.ty, 2 * design.taupm;
                          -outermagsep/2 + design.g + design.tc, 2 * design.taupm ];

            % add the nodes at the coordinates
            [FemmProblem, yokenodeinds, yokenodeids] = addnodes_mfemm(FemmProblem, yokecoords(:,1), yokecoords(:,2));

            % get the location for a yoke label for later
            yokelabelloc = rectcentre(yokecoords(1,:), yokecoords(3,:));
            
            % add a circuit for the coils
            FemmProblem = addcircuit_mfemm(FemmProblem, 'A', ...
                'TotalAmps_re', real(Inputs.CoilCurrent), ...
                'TotalAmps_im', imag(Inputs.CoilCurrent));

            % add the coil material which should be present in the design
            % structure
            FemmProblem.Materials = [FemmProblem.Materials, matstr2matstruct_mfemm(design.MagFEASimMaterials.ArmatureCoil)];

            % define the block properties of the coil region
            coilBlockProps.BlockType = FemmProblem.Materials(end).Name;
            coilBlockProps.MaxArea = Inputs.AirGapMeshSize;
            coilBlockProps.InCircuit = 'A';
            coilBlockProps.InGroup = Inputs.CoilGroup;

            % draw the positive part of the coil circuit
            coilBlockProps.Turns = design.CoilTurns;

            xcoil = -outermagsep/2 + design.g;
            ycoil = design.taupm - design.tauco/2;
            
            % add a rectangular region making up half the coil
            [FemmProblem, coilseginds, coilnodeinds, coilblockind, coilnodeids] = ...
                addrectregion_mfemm(FemmProblem, xcoil, ycoil, design.tc, design.tauco, coilBlockProps); 

            % link up the top and bottom of the yoke at the coil edges
            
            % the top of coil to top of yoke
            FemmProblem = addsegments_mfemm(FemmProblem, coilnodeids(3), yokenodeids(4));
            
            % the bottom of coil to bottom of yoke
            FemmProblem = addsegments_mfemm(FemmProblem, coilnodeids(2), yokenodeids(1));
            
            % draw the -ve part of the coil circuit
            coilBlockProps.Turns = -design.CoilTurns;

            xcoil = xcoil + design.tc + design.ty;

            % add a rectangular region making up the other half of the coil
            [FemmProblem, coilseginds, coilnodeinds, coilblockind, coilnodeids] = ...
                addrectregion_mfemm(FemmProblem, xcoil, ycoil, design.tc, design.tauco, coilBlockProps);        

            % link up the top and bottom of the yoke at the coil edges

            % the bottom of coil to bottom of yoke
            FemmProblem = addsegments_mfemm(FemmProblem, coilnodeids(1), yokenodeids(2));
            
            % the top of coil to top of yoke
            FemmProblem = addsegments_mfemm(FemmProblem, coilnodeids(4), yokenodeids(3));
            
            % Add segments to the top and bottom, with periodic boundary
            % conditions
            
            % the top of the yoke
            FemmProblem = addsegments_mfemm(FemmProblem, yokenodeids(1), yokenodeids(2), SegProps(1));
            
            % the bottom
            FemmProblem = addsegments_mfemm(FemmProblem, yokenodeids(4), yokenodeids(3), SegProps(3));
            
            % add a yoke label
            FemmProblem = addblocklabel_mfemm(FemmProblem, yokelabelloc(1,1), yokelabelloc(1,2), struct2pvpairs(yokeBlockProps));
            
        else

            % add the core region
            [FemmProblem, seginds, yokenodeinds, blockind, yokenodeids] = ...
                addrectregion_mfemm(FemmProblem, corexpos, 0, design.ty, 2 * design.taupm, yokeBlockProps, SegProps);
        
        end
        
        % join up the air gaps at the top and bottom
        
        % bottom left gap corner to bottom left core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(1), yokenodeids(1), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(1)).Name);

        % top left gap corner to top left core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(3), yokenodeids(4), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(1)).Name);   

        % bottom right gap corner to bottom right core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(2), yokenodeids(2), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);

        % top right gap corner to top right core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(4), yokenodeids(3), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);  

        corexpos = corexpos + innerstagewidth;
        
        % Add block labels for the air gaps
        labelloc = [gapedgenodes(1,1) + design.g/2, design.taupm];

        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', FemmProblem.Materials(1).Name, ...
                                        'MaxArea', Inputs.AirGapMeshSize);
                                    
        labelloc = [gapedgenodes(1,1) + 3*design.g/2 + 2*design.tc + design.ty, design.taupm];

        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', FemmProblem.Materials(1).Name, ...
                                        'MaxArea', Inputs.AirGapMeshSize);
        
        % shift the nodes to the next location
        gapedgenodes(:,1) = gapedgenodes(:,1) + innerstagewidth;
        
    end

end

% 
% function design = checkdesign_()
% 
% SegProps  valfields = {'Rrmo', 'Rrmi', 'taupro', 'tauprm', 'Rsci', 'Rsco', 'taupso', 'tTaupsi'};
%     
%     ratiofields = {};
%     
%     basefield = {'Ro', };
%     
%     design = checkdesign_am(design, valfields, ratiofields, basefield, varargin)
% 
% end
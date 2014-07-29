function [FemmProblem, outermagsep] = corelessLfemmprob_torus(design, varargin)
% corelesstorusfemmprob creates a FemmProblem structure for a coreless (or
% yokeless) torus axial flux permanent magnet machine

    magsep = 2*design.g + design.tc;
    
    Inputs.NStages = 1;
    Inputs.DrawCoils = true;
    Inputs.CoilCurrent = 0;
    Inputs.MagArrangement = 'NS';
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.Rmo - design.Rmi);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
    Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
    Inputs.MagnetSpaceGroup =[];
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
    
    FemmProblem = Inputs.FemmProblem;
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    FemmProblem = addmaterials_mfemm( FemmProblem, ...
                             matstr2matstruct_mfemm( {design.MagSimMaterials.Magnet, ...
                                                      design.MagSimMaterials.FieldIron} ));
    
    % draw the torus rotor according to the spec in the design strucure
    [FemmProblem, outermagsep, innerstagewidth] = torusrotor2dfemmprob( ...
        design.taupm, design.taumm, design.tm, design.tbi, magsep, ...
        'NStages', Inputs.NStages, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', 1, ... % replace mags with air
        'BackIronMaterial', numel(FemmProblem.Materials), ...
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
    
    % link the rotor stages along the top and bottom, and add antiperiodic
    % boundaries

    
    gapedgenodes = [-outermagsep/2, 0;
                    -outermagsep/2+magsep, 0;
                    -outermagsep/2, 2*design.taupm;
                    -outermagsep/2+magsep, 2*design.taupm];
    
    for i = 1:Inputs.NStages
        
        % get the node ids of the air gap corners
        nodeids = findnode_mfemm(FemmProblem, gapedgenodes);
        
        % add a new periodic boundary for the top and bottom of the region
        [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Rect Mags Air Gap Periodic', 4);

        FemmProblem.Segments(end+1) = newsegment_mfemm(nodeids(1), nodeids(2), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name);

        FemmProblem.Segments(end+1) = newsegment_mfemm(nodeids(3), nodeids(4), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind).Name);   

        % Add block label for the back iron
        labelloc = [gapedgenodes(1,1) + design.g/2, design.taupm];

        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', FemmProblem.Materials(1).Name, ...
                                        'MaxArea', Inputs.AirGapMeshSize);
        
        % shift the nodes to the next location
        gapedgenodes(:,1) = gapedgenodes(:,1) + innerstagewidth;
        
    end
    
    % now modify the drawing to encompass a set of adjacent Phases
    
    if strcmp(design.WindingType, 'nonoverlapping')
        % scale all the node and block coordingate in the vertical direction to
        % give the sim height equal to a repeating set of adjacent Phases
        FemmProblem = scaleproblem_mfemm(FemmProblem, 1, design.taupcg / (2*design.taupm));
        for i = 1:numel(FemmProblem.BlockLabels)
            FemmProblem.BlockLabels(i).MaxArea = FemmProblem.BlockLabels(i).MaxArea * design.taupcg / (2*design.taupm);
        end
        
        ycoil = design.taupcg/design.Phases/2;
        
        ycoildisp = design.taupcg/design.Phases;
        
    elseif strcmp(design.WindingType, 'overlapping')
        
%         FemmProblem = scaleproblem_mfemm(FemmProblem, 1, design.taupm / (2*design.taupm));
        
%         for i = 1:numel(FemmProblem.BlockLabels)
%             FemmProblem.BlockLabels(i).MaxArea = FemmProblem.BlockLabels(i).MaxArea * design.taupm / (2*design.taupm);
%         end
        
        ycoil = 3*design.taupm/4;
        
        ycoildisp = design.taupm/design.Phases;
        
    else
        error('CORELESS_TORUS:corelessLfemmprob:badwindingtype', ...
              'Unrecognised winding type.')
    end
    
    % Now draw the coils 
    for i = 1:2 %design.Phases
        
        % add a circuit for the coil
        if i == 1
            % only apply current to the first coil
            FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i), ...
                'TotalAmps_re', real(Inputs.CoilCurrent), ...
                'TotalAmps_im', imag(Inputs.CoilCurrent));
        else
            FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i), ...
                'TotalAmps_re', 0, ...
                'TotalAmps_im', 0);
        end
        
        % add the coil material which should be present in the design
        % structure
        FemmProblem.Materials = [FemmProblem.Materials, matstr2matstruct_mfemm(design.MagSimMaterials.CoilWinding)];
        
        % define the block properties of the coil region
        BlockProps.BlockType = FemmProblem.Materials(end).Name;
        BlockProps.MaxArea = Inputs.AirGapMeshSize;
        BlockProps.InCircuit = num2str(i);
        BlockProps.InGroup = Inputs.CoilGroup;
        
        % draw the positive part of the coil circuit
        BlockProps.Turns = design.CoilTurns;
    
        xcoil = -outermagsep/2 + design.g;
        ycoiltop = ycoil + design.tauci/2;
        
        % add a rectangular region making up half the coil
        FemmProblem = addrectregion_mfemm(FemmProblem, xcoil, ycoiltop, ...
            design.tc/design.CoilLayers, (design.tauco - design.tauci) / 2, BlockProps);
        
        % draw the -ve part of the coil circuit
        BlockProps.Turns = -design.CoilTurns;
    
        xcoil = xcoil + (design.tc/design.CoilLayers)*(design.CoilLayers-1);
        ycoilbot = ycoil - design.tauco/2;
        
        % add a rectangular region making up the other half of the coil
        FemmProblem = addrectregion_mfemm(FemmProblem, xcoil, ycoilbot, ...
            design.tc/design.CoilLayers, (design.tauco - design.tauci) / 2, BlockProps);  
        
        ycoil = ycoil + ycoildisp;
        
    end

end

% 
% function design = checkdesign_()
% 
%     valfields = {'Rrmo', 'Rrmi', 'taupro', 'tauprm', 'Rsci', 'Rsco', 'taupso', 'taupsi'};
%     
%     ratiofields = {};
%     
%     basefield = {'Ro', };
%     
%     design = checkdesign_am(design, valfields, ratiofields, basefield, varargin)
% 
% end
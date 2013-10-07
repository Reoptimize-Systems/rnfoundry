function [FemmProblem, outermagsep] = slotlesseddyfemmprob_torus(design, RPM, varargin)
% creates a FemmProblem structure for a slotless torus axial flux permanent
% magnet machine

    magsep = 2*design.g + 2*design.tc + design.ty;
    
    Inputs.NStages = 1;
    Inputs.CoilCurrent = 0;
    Inputs.MagArrangement = 'NN';
    Inputs.FemmProblem = newproblem_mfemm('planar', 'Depth', design.Rmo - design.Rmi, 'freq',  RPM / (60 * design.Poles) );
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
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    elcount = elementcount_mfemm(FemmProblem);
    
    % Convert the material names to materials structures from the materials
    % library, if this has not already been done.
    FemmProblem.Materials = [FemmProblem.Materials, ...
                             matstr2matstruct_mfemm( {design.MagnetMaterial, ...
                                                      design.BackIronMaterial} )];
    
    BackironMatInd = elcount.NMaterials + 2;
    
    % draw the torus rotor according to the spec in the design strucure
    [FemmProblem, outermagsep, innerstagewidth] = torusrotor2dfemmprob( ...
        design.taupm, design.taumm, design.tm, design.tbi, magsep, ...
        'NStages', Inputs.NStages, ...
        'FemmProblem', FemmProblem, ...
        'MagArrangement', Inputs.MagArrangement, ...
        'MagnetMaterial', 1, ...
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
    magnxpos = -outermagsep/2 - design.tm;
    
    ycurrentpoints = 10;
    xcurrentsecs = 5;
    
    yn = linspace(design.taupm / (2*ycurrentpoints), (2 * design.taupm) - (design.taupm / (2*ycurrentpoints)), ycurrentpoints);
    
    xn = linspace((design.tm / xcurrentsecs)/2, design.tm - ((design.tm / xcurrentsecs)/2), xcurrentsecs);
    
    J = design.HcMag / (design.taupm / ycurrentpoints);
     
    pointcurrentmag = J * (design.taupm / ycurrentpoints) * (design.tm / xcurrentsecs);
    
    pointcurrents = pointcurrentmag * (cos(pi * yn./design.taupm) + 1j*sin(pi * yn./design.taupm));
    
    % add the current point properties for later use
    for i = 1:numel(yn)
        
        PointProp.I_re = real(pointcurrents(i));
        PointProp.I_im = imag(pointcurrents(i));
        
        [FemmProblem, ppropinds(i)] = addpointprop_mfemm(FemmProblem, int2str(i), PointProp);
        
    end
    
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
        
        % add the core region
        [FemmProblem, seginds, yokenodeinds, blockind, yokenodeids] = ...
            addrectregion_mfemm(FemmProblem, corexpos, 0, design.ty, 2 * design.taupm, yokeBlockProps, SegProps);
        
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
        
        % now add the travelling wave of currents representing the
        % permanent magnets
        
        % add left hand magnets
        for ii = 1:numel(yn)
            
            FemmProblem = addnodes_mfemm(FemmProblem, xn(:)+magnxpos, repmat(yn(ii), numel(xn), 1), 'PointPropName', int2str(ii));
            
        end
        
        magmove = design.tm + design.tm + 2*(design.g + design.tc);
        
        % add right hand magnets
        switch Inputs.MagArrangement
            
            case 'NN'
                
                for ii = 1:numel(yn)
                    
                    FemmProblem = addnodes_mfemm(FemmProblem, xn(:)+magnxpos+magmove, repmat(yn(ii), numel(xn), 1), 'PointPropName', int2str(numel(yn)+1-ii));
                    
                end
                
            case 'NS'
                
                for ii = 1:numel(yn)
                    
                    FemmProblem = addnodes_mfemm(FemmProblem, xn(:)+magnxpos+magmove, repmat(yn(ii), numel(xn), 1), 'PointPropName', int2str(ii));
                    
                end
                
            otherwise
                
        end
        
        % shift the nodes to the next location
        gapedgenodes(:,1) = gapedgenodes(:,1) + innerstagewidth;
        magnxpos = magnxpos + innerstagewidth;
        
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
function [FemmProblem, fieldinfo, statorinfo] = slotlessfemmprob_tubular(design, varargin)
% creates a FemmProblem structure for a slotless radial flux permanent
% magnet machine
%
% Syntax
%
% [FemmProblem, rotorinfo, statorinfo] = slotlessfemmprob_radial(design)
% [FemmProblem, rotorinfo, statorinfo] = slotlessfemmprob_radial(design, 'Parameter', Value)
%
% 


    % set some parameters which may not be present in the design structure
    design = setfieldifabsent (design, 'rsb', 0);
    design = setfieldifabsent (design, 'rsg', 0);
    design = setfieldifabsent (design, 'zsg', 0);
               
    % *****                 Function input parsing                 ****** %
%     Inputs.DrawingType = 'MagnetRotation';
    Inputs.NBoundaryPositions = 10;
    Inputs.BoundaryShift = 0;
    Inputs.NWindingLayers = nan;
    Inputs.CoilCurrent = zeros (1,design.Phases);
    Inputs.MagArrangement = 'NN';
    Inputs.PolarisationType = 'yz';
    if isfield (design, 'MagnetPolarisation') && ischar (design.MagnetPolarisation)
        Inputs.PolarisationType = design.MagnetPolarisation;
    end
    Inputs.FemmProblem = newproblem_mfemm ('axi', 'MinAngle', 15);
    Inputs.Position = 0;
    Inputs.FractionalPolePosition = [];
%     Inputs.RotorAnglePosition = [];
    Inputs.MagnetGroup = [];
%     Inputs.MagnetSpaceGroup = [];
%     Inputs.FieldBackIronGroup = [];
    Inputs.FieldOuterRegionGroup = [];
    Inputs.CoilGroup = 0;
    Inputs.ArmatureBackIronGroup = [];
    Inputs.MagnetRegionMeshSize = choosemesharea_mfemm (design.rm, design.zm, 1/10);
    Inputs.MagnetSpacerRegionMeshSize = choosemesharea_mfemm (design.rm, design.zs, 1/10);
%     Inputs.FieldIronRegionMeshSize = choosemesharea_mfemm (design.rbi, 2*design.zp, 1/10);
    Inputs.FieldOuterRegionsMeshSize = [choosemesharea_mfemm(design.rm, design.zp, 1/5), -1];
    Inputs.StatorOuterRegionsMeshSize = [];
    Inputs.StatorOuterRegionMaterials = {};
    Inputs.AirGapMeshSize = choosemesharea_mfemm (design.g, design.zp, 1/10);
    Inputs.DrawOuterRegions = true;
    Inputs.StatorOuterRegionSize = [];
    Inputs.CoilInsRegionMeshSize = -1;
    
    if design.rsg > 1e-5
        if design.rsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (max(design.rsg, design.rsb), design.zsg, 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (design.rsb, design.zsg, 1/20);
        end
    else
        if design.rsb > 1e-5
            Inputs.ShoeGapRegionMeshSize = ...
                choosemesharea_mfemm (design.rsb, design.zsg, 1/20);
        else
            Inputs.ShoeGapRegionMeshSize = -1;
        end
    end
    Inputs.YokeRegionMeshSize = mean( [choosemesharea_mfemm(design.ry, 2*design.zp, 1/10), ...
                                       choosemesharea_mfemm(design.rc(1), design.zs-mean(design.zc), 1/10)] );
    Inputs.CoilRegionMeshSize = choosemesharea_mfemm (design.rc(1), mean(design.zc));
    Inputs.Tol = 1e-5;
    Inputs.SimType = 'Magnetics';
    Inputs.MaterialsLibrary = '';
    Inputs.NPolePairs = 1;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    % never draw coil insulation
    Inputs.DrawCoilInsulation = false;
    
    % *****                 Function input parsing                 ****** %
    
    % we'll draw a radial flux stator with the stator iron replaced with
    % air
    
    switch design.ArmatureType
        
        case 'external'
            % single inner facing stator
            Rs = design.Rag + design.g + design.rc(1) + design.rsb + design.ry/2;
            
            if isempty (Inputs.StatorOuterRegionSize)
                Inputs.StatorOuterRegionSize = [2*design.rm, 10*design.rm];
            end
            
        case 'internal'
            
            % single outer facing stator
            Rs = design.Rag - design.g - design.rc(1) - design.rsb - design.ry/2;
            
            if isempty (Inputs.StatorOuterRegionSize)
                Inputs.StatorOuterRegionSize = [0.8, 0.5];  
            end
            
        otherwise
            error('Unrecognised ArmatureType option.')
                
    end
    
    if isempty (Inputs.StatorOuterRegionsMeshSize)
        Inputs.StatorOuterRegionsMeshSize = [ choosemesharea_mfemm(design.rm, design.zp, 1/5), ...
                                              repmat(-1,1,numel(Inputs.StatorOuterRegionSize)-1) ];
    end
    
    if isempty (Inputs.StatorOuterRegionMaterials)
        Inputs.StatorOuterRegionMaterials = ...
                repmat({design.MagFEASimMaterials.AirGap}, 1, numel(Inputs.StatorOuterRegionSize));
    end
    
    % do some checking 
    assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionMaterials), ...
        'RENEWNET:slottedfemmproblem_radial:nstatormaterials', ...
        'Number of supplied stator outer region material names does not match number of outer region sizes.');
    
    assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionsMeshSize), ...
                'RENEWNET:slotlessfemmproblem_radial:nstatormeshsizes', ...
                'Number of supplied stator outer region mesh sizes does not match number of outer region sizes.');
    
    if (numel (design.rc) == 1) || (design.rc(2) <= Inputs.Tol)
        % make it zero
        design.rc(2) = 0;
        
    else
        
        % make the yoke drawing thickness as small as possible, we will be
        % using the first outer region as the yoke
        design.ty = min (1e-3, 5 * Inputs.Tol);
        
        switch design.ArmatureType
            
            case 'external'
                Inputs.StatorOuterRegionSize = [design.ty, Inputs.StatorOuterRegionSize];
                
            case 'internal'
                Inputs.StatorOuterRegionSize = [Inputs.StatorOuterRegionSize, design.Ryi / design.Ryo];
                
        end
        
        Inputs.StatorOuterRegionsMeshSize = [Inputs.YokeRegionMeshSize, Inputs.StatorOuterRegionsMeshSize];        
        
        % handle the outer region materials
        % get what the armature yoke material should be and store it
        yokemagsimmat = design.MagFEASimMaterials.ArmatureYoke;
        % swap the yoke material to be the same as the air gap material
        design.MagFEASimMaterials.ArmatureYoke = design.MagFEASimMaterials.AirGap;
        
        Inputs.StatorOuterRegionMaterials = [ yokemagsimmat, Inputs.StatorOuterRegionMaterials ];
    
    end
    
    assert (samesize (Inputs.StatorOuterRegionSize, Inputs.StatorOuterRegionMaterials), ...
        'RENEWNET:slotlessfemmproble_radial:nstatormaterials', ...
        'Number of supplied stator outer region material names does not match number of outer region sizes.');
        
    % draw the iron cored machine
    [FemmProblem, fieldinfo, statorinfo] = slottedfemmprob_tubular ( design, Inputs );
    
    if design.rc(2) == 0
        % there is no curved slot/coil base, so we need to link up the
        % corners of adjacent slots to form the inner boundary of the yoke
       
        % calculate slot corner node locations
        ycLower = linspace (-design.zcy/2, ...
                            (statorinfo.NDrawnSlots-1) * design.zs - design.zcy/2, ...
                                statorinfo.NDrawnSlots ) ...
                          + design.zs/2;
                      
        ycUpper = ycLower + design.zcy;
        
        ycLower = ycLower + statorinfo.ZShift;
        ycUpper = ycUpper + statorinfo.ZShift;
        
        xcLower = repmat (design.Ryi, size (ycLower));
        xcUpper = repmat (design.Ryi, size (ycUpper));
        
        % link up slot sides using arc segments
        nidLower = findnode_mfemm (FemmProblem, [xcLower(:), ycLower(:)]);
        nidUpper = findnode_mfemm (FemmProblem, [xcUpper(:), ycUpper(:)]);
        
        zToothWidth = design.zs - design.zcy;
        
        % first the slots to each other
        FemmProblem = addsegments_mfemm (FemmProblem, ...
                                         nidUpper(1:end-1), nidLower(2:end));
        
        % add block labels for teeth
        ytlabel = ycUpper + zToothWidth/2;
                            
        xtlabel = repmat (design.Rcm, 1, numel (ytlabel)-1);
                            
        for ind = 1:numel (xtlabel)
            % we add all but the last tooth label
            FemmProblem = addblocklabel_mfemm (FemmProblem, xtlabel(ind), ytlabel(ind), ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );
        end
        
        % next handle the bottom and top slots, what is done here depends
        % on whether the full machine is being drawn or not
        if statorinfo.NDrawnSlots == design.Qs
            % join the last slot to the first slot
            FemmProblem = addsegments_mfemm ( FemmProblem, ...
                                              nidUpper(end), nidLower(1) );
            
            % add a single label for the last tooth
            FemmProblem = addblocklabel_mfemm ( FemmProblem, design.Rcm, statorinfo.ZShift, ...
                                                'BlockType', design.MagFEASimMaterials.AirGap );
            
        else
            % join the top slot to the top and the bottom slot to the
            % bottom
            
            % find and the split the top and bottom segments to create the
            % required nodes at the appropriate positions
            segx = [design.Ryi; design.Ryi];
            segy = [0; design.zp*Inputs.NPolePairs*2] + statorinfo.ZShift;
            segids = findsegment_mfemm (FemmProblem, [segx, segy]);
            
            seglength = segmentlength_mfemm (FemmProblem, segids(1));
            
            segstartcoords = FemmProblem.Nodes(FemmProblem.Segments(segids(1)+1).n0+1).Coords;
            segendcoords = FemmProblem.Nodes(FemmProblem.Segments(segids(1)+1).n1+1).Coords;
            
            lenfrac = ( design.Ryi - min (segstartcoords(1), segendcoords(1)) ) / seglength;
            
            [FemmProblem, newbotsegid, newbotnodeid] = splitsegment_mfemm (FemmProblem, segids(1), lenfrac);
            [FemmProblem, newtopsegid, newtopnodeid] = splitsegment_mfemm (FemmProblem, segids(2), lenfrac);
            
            % add new periodic boundaries top and bottom to the newly
            % created segments
            [FemmProblem, ~, boundname] = addboundaryprop_mfemm (FemmProblem, ...
                                                                 'Periodic', 4);
                                                             
            % swap the boundary marker to avoid confusion (it's name will
            % say it's the back iron boundary marker)
            FemmProblem.Segments(newbotsegid+1).BoundaryMarker = FemmProblem.Segments(segids(1)+1).BoundaryMarker;
            FemmProblem.Segments(newtopsegid+1).BoundaryMarker = FemmProblem.Segments(segids(2)+1).BoundaryMarker;
            % add the new boundary marker
            FemmProblem.Segments(segids(1)+1).BoundaryMarker = boundname;
            FemmProblem.Segments(segids(2)+1).BoundaryMarker = boundname;
            
            % join the bottom slot to the new bottom node
            FemmProblem = addsegments_mfemm (FemmProblem, ...
                                             newbotnodeid, nidLower(1));
            
            % join the top slot to the new top node
            FemmProblem = addsegments_mfemm (FemmProblem, ...
                                             nidUpper(end), newtopnodeid );

            % add two labels, one for each tooth half
            x = design.Rcm;
            y = (design.zs-max(design.zcy, design.zcg))/4 ...
                + statorinfo.ZShift;
            
            FemmProblem = addblocklabel_mfemm (FemmProblem, x, y, ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );

            x = design.Rcm;
            y = design.zp*Inputs.NPolePairs*2 - (design.zs-max(design.zcy, design.zcg))/4 ...
                + statorinfo.ZShift;
            
            FemmProblem = addblocklabel_mfemm (FemmProblem, x, y, ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );
            
        end
        
    end
    
end


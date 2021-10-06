function [FemmProblem, rotorinfo, statorinfo] = slotlessfemmprob_radial(design, varargin)
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
    design = setfieldifabsent (design, 'tsb', 0);
    design = setfieldifabsent (design, 'tsg', 0);
    design = setfieldifabsent (design, 'thetasg', 0);
               
    % *****                 Function input parsing                 ****** %
    Inputs.DrawingType = 'MagnetRotation';
    Inputs.NBoundaryPositions = 10;
    Inputs.BoundaryShift = 0;
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

    Inputs = parse_pv_pairs (Inputs, varargin);
    
    % never draw coil insulation
    Inputs.DrawCoilInsulation = false;
    
    % *****                 Function input parsing                 ****** %
    
    % we'll draw a radial flux stator with the stator iron replaced with
    % air
    
    switch design.ArmatureType
        
        case 'external'
            % single inner facing stator
            Rs = design.Rmo + design.g + design.tc(1) + design.tsb + design.ty/2;
            
            if isempty (Inputs.StatorOuterRegionSize)
                Inputs.StatorOuterRegionSize = [2*design.tm, 10*design.tm];
            end
            
        case 'internal'
            
            % single outer facing stator
            Rs = design.Rmi - design.g - design.tc(1) - design.tsb - design.ty/2;
            
            if isempty (Inputs.StatorOuterRegionSize)
                Inputs.StatorOuterRegionSize = [0.8, 0.5];  
            end
            
            
        case 'di'
            % double internal stator (mags on outside)
%             Rs = design.Rmo(1) + design.g + design.tc(1) + design.ty/2;
            error('not yet supported');
        case 'do'
            % double outer/external stator (mags on inside)
            error('not yet supported');
            
        otherwise
            error('Unrecognised ArmatureType option.')
                
    end
    
    if isempty (Inputs.StatorOuterRegionsMeshSize)
        Inputs.StatorOuterRegionsMeshSize = [ choosemesharea_mfemm(design.tm, (Rs*design.thetap), 1/5), ...
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
    
    if (numel (design.tc) == 1) || (design.tc(2) <= Inputs.Tol)
        % make it zero
        design.tc(2) = 0;
        
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
    [FemmProblem, rotorinfo, statorinfo] = slottedfemmprob_radial ( design, Inputs );
    
    if design.tc(2) == 0
        % there is no curved slot/coil base, so we need to link up the
        % corners of adjacent slots to form the inner boundary of the yoke
       
        % calculate slot corner node locations
        thetaLower = linspace (-design.thetacy/2, ...
                                (statorinfo.NDrawnSlots-1)* design.thetas - design.thetacy/2, ...
                                statorinfo.NDrawnSlots ) + design.thetas/2;
        thetaUpper = thetaLower + design.thetacy;
        
        [xcLower, ycLower] = pol2cart (thetaLower, repmat (design.Ryi, size (thetaLower)));
        [xcUpper, ycUpper] = pol2cart (thetaUpper, repmat (design.Ryi, size (thetaUpper)));
        
        % link up slot sides using arc segments
        nidLower = findnode_mfemm (FemmProblem, [xcLower(:), ycLower(:)]);
        nidUpper = findnode_mfemm (FemmProblem, [xcUpper(:), ycUpper(:)]);
        
        thetaToothWidth = design.thetas - design.thetacy;
        
        % first the slots to each other
        FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                            nidUpper(1:end-1), nidLower(2:end), ...
                                            repmat (rad2deg (thetaToothWidth), 1, numel (nidUpper)-1));
        
        % add block labels for teeth
        thetaToothMiddle = thetaUpper + thetaToothWidth/2;
                            
        [xtlabel, ytlabel] = pol2cart (thetaToothMiddle(1:end-1), ...
                                       repmat (design.Rcm, 1, numel (thetaToothMiddle)-1));
                            
        for ind = 1:numel (xtlabel)
            % we add all but the last tooth label
            FemmProblem = addblocklabel_mfemm (FemmProblem, xtlabel(ind), ytlabel(ind), ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );
        end
        
        % next handle the bottom and top slots, what is done here depends
        % on whether the full machine is being drawn or not
        if statorinfo.NDrawnSlots == design.Qs
            % join the last slot to the first slot
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                nidUpper(end), nidLower(1), ...
                                                rad2deg (thetaToothWidth) );
            
            % add a single label for the last tooth
            FemmProblem = addblocklabel_mfemm (FemmProblem, design.Rcm, 0, ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );
            
        else
            % join the top slot to the top and the bottom slot to the
            % bottom
            
            % find and the split the top and bottom segments to create the
            % required nodes at the appropriate positions
            [segx, segy] = pol2cart ([0; design.thetap*Inputs.NPolePairs*2], [design.Ryi; design.Ryi]);
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
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                newbotnodeid, nidLower(1), ...
                                                rad2deg (thetaToothWidth)/2 );
            
            % join the top slot to the new top node
            FemmProblem = addarcsegments_mfemm (FemmProblem, ...
                                                nidUpper(end), newtopnodeid, ...
                                                rad2deg (thetaToothWidth)/2 );

            % add two labels, one for each tooth half
            
            [x, y] = pol2cart ((design.thetas-max(design.thetacy, design.thetacg))/4, design.Rcm);
            
            FemmProblem = addblocklabel_mfemm (FemmProblem, x, y, ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );

            [x, y] = pol2cart (design.thetap*Inputs.NPolePairs*2 - (design.thetas-max(design.thetacy, design.thetacg))/4, design.Rcm);
            
            FemmProblem = addblocklabel_mfemm (FemmProblem, x, y, ...
                                               'BlockType', design.MagFEASimMaterials.AirGap );
            
        end
        
    end
    
end


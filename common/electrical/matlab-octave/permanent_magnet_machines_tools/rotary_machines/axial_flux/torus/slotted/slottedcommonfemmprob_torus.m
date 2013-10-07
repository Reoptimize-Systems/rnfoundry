function FemmProblem = slottedcommonfemmprob_torus(FemmProblem, design, ...
    Inputs, magsep, gapedgenodes, innerstagewidth, coillabellocs, yokenodeids, ...
    ymidpoint, BackIronMatInd, YokeMatInd, CoilMatInd)

    slotsperpole = design.Qs / design.Poles;

    % define the block properties of the core region
    yokeBlockProps.BlockType = FemmProblem.Materials(BackIronMatInd).Name;
    yokeBlockProps.MaxArea = Inputs.BackIronRegionMeshSize;
    yokeBlockProps.InCircuit = '';
    yokeBlockProps.InGroup = 0;
        
    % Prototype an array of segprops structures
    SegProps.MaxSideLength = -1;
    SegProps.Hidden = 0;
    SegProps.InGroup = 0;
    SegProps.BoundaryMarker = '';
    
    SegProps = repmat(SegProps, 1, 4);
    
    coilBlockProps.BlockType = FemmProblem.Materials(CoilMatInd).Name;
    coilBlockProps.MaxArea = Inputs.CoilRegionMeshSize;
    coilBlockProps.InCircuit = '';
    coilBlockProps.InGroup = Inputs.CoilGroup;

    % draw the positive part of the coil circuit
    coilBlockProps.Turns = design.CoilTurns;
    
%     corexpos = -outermagsep/2 + design.g + design.tc;

    % add circuits for each phase
    for i = 1:design.Phases
        if i == 1
            FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i), 'TotalAmps_re', Inputs.CoilCurrent);
        else
            FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i));
        end
    end
    
    for i = 1:Inputs.NStages
        
        % get the node ids of the air gap corners on the rotor stages
        gapcornernodeids = findnode_mfemm(FemmProblem, gapedgenodes);
        
        % add four new nodes at the boundary of the air gap and the teeth
        % on the armature
        newgapnodes = [gapedgenodes(1,1) + design.g, gapedgenodes(1,2);
                       gapedgenodes(2,1) - design.g, gapedgenodes(2,2);
                       gapedgenodes(3,1) - design.g, gapedgenodes(3,2)
                       gapedgenodes(4,1) + design.g, gapedgenodes(4,2)];
                   
        [FemmProblem, outeryokenodeinds, outeryokenodeids] = ...
            addnodes_mfemm(FemmProblem, newgapnodes(:,1), newgapnodes(:,2));
        
        % add a three new periodic boundaries for the top and bottom of the
        % air gap and core regions
        [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Stator Air Gap Periodic', 4);
        [FemmProblem, boundind(2)] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Stator Back Iron Periodic', 4);
        [FemmProblem, boundind(3)] = addboundaryprop_mfemm(FemmProblem, 'Multi Stage Stator Air Gap Periodic', 4);
        
        % add the segments for the top and bottom of the core       
        [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, ...
                                    [outeryokenodeids(1); outeryokenodeids(4)], ...
                                    [outeryokenodeids(2); outeryokenodeids(3)], ...
                                    'BoundaryMarker', FemmProblem.BoundaryProps(boundind(2)).Name);
                
        % join up the air gaps at the top and bottom
        
        % bottom left gap corner to bottom left core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(1), outeryokenodeids(1), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(1)).Name);

        % top left gap corner to top left core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(4), outeryokenodeids(4), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(1)).Name);   

        % bottom right gap corner to bottom right core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(2), outeryokenodeids(2), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);

        % top right gap corner to top right core corner
        FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(3), outeryokenodeids(3), ...
                                'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);  

        % close the tooth and air boundaries to complete the core
        FemmProblem = addsegments_mfemm(FemmProblem, yokenodeids(i,1:4), outeryokenodeids(1:4));  
                            
%         corexpos = corexpos + innerstagewidth;
        
        % Add block labels for the air gaps
        labelloc = [gapedgenodes(1,1) + design.g/2, ymidpoint];

        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', FemmProblem.Materials(1).Name, ...
                                        'MaxArea', Inputs.AirGapMeshSize);
                                    
        labelloc = [gapedgenodes(2,1) - design.g/2, ymidpoint];

        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', FemmProblem.Materials(1).Name, ...
                                        'MaxArea', Inputs.AirGapMeshSize);
                                    
        % add a block label for the yoke and teeth
        labelloc = [gapedgenodes(1,1) + magsep/2, ymidpoint];

        FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                        'BlockType', FemmProblem.Materials(YokeMatInd).Name, ...
                                        'MaxArea', Inputs.YokeRegionMeshSize);
        
        % add block labels for the coils
        j = 2*(i-1) + 1;

        row = 1;
        
        circnums = zeros(Inputs.NSlots, 1);
        temp = (1:design.Phases)';
        
        if design.yd == 1
            
            % short pitched winding
            for ii = 1:2:Inputs.NSlots

                circnums(ii) = temp(1);
                
                if  ii < numel(circnums)
                    
                    circnums(ii+1) = temp(1);
                    
                end

                temp = circshift(temp, 1);
                
            end
        
        else
            
            % otherwise next slot contains the next phase, and so on in
            % sequence
            for ii = 1:Inputs.NSlots

                circnums(ii) = temp(1);

                temp = circshift(temp, 1);
                
            end
            
            
        end
        
        docircname = zeros(numel(circnums), Inputs.NWindingLayers);
        
        if design.yd == 1
            
            for ii = 1:2:2*design.Phases

                if ii <= numel(circnums)

                    docircname(ii, :) = [1, zeros(1, Inputs.NWindingLayers-1)];

                end

                if ii+design.yd <= numel(circnums)

                    docircname(ii+design.yd, :) = [zeros(1, Inputs.NWindingLayers-1), -1];

                end

            end

        else

            for ii = 1:design.Phases

                if ii <= numel(circnums)

                    docircname(ii, :) = [1, zeros(1, Inputs.NWindingLayers-1)];

                end

                if ii+design.yd <= numel(circnums)

                    docircname(ii+design.yd, :) = [zeros(1, Inputs.NWindingLayers-1), -1];

                end

            end

        end

%         circslotcount = 1;
        
%         slotnums = (1:Inputs.NSlots)';
%         nextlayer = 1;
%         layers = (1:Inputs.NWindingLayers)';
        
        for k = 1:Inputs.NSlots

            for n = 1:Inputs.NWindingLayers

                if k <= 2*design.Phases && docircname(k,n) ~= 0
                    
                    % only put the circuit in the first set of phase coils
                    coilBlockProps.InCircuit = num2str(circnums(k));
                    coilBlockProps.Turns = coilBlockProps.Turns * docircname(k,n);
                
                else
                    
                    % only put the circuit in the first set of phase coils
                    coilBlockProps.InCircuit = '';
                    
                end 
                
                FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                    coillabellocs(row,j), coillabellocs(row,j+1), ...
                    coilBlockProps);

%                 row = row + 1;

                FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                    coillabellocs(row+(Inputs.NSlots*Inputs.NWindingLayers),j), coillabellocs(row,j+1), ...
                    coilBlockProps);
                
                coilBlockProps.Turns = abs(coilBlockProps.Turns);
                coilBlockProps.InCircuit = '';

                row = row + 1;

            end
            
%             nextlayer = circshift(layers, -1);
                
%             nextlayer = nextlayer(1);
            
%             circslotcount = circslotcount + 1;
            
%             circnums = circshift(circnums, 1);
            
%             slotnums = circshift(slotnums, );

        end
        
        % shift the nodes to the next location
        gapedgenodes(:,1) = gapedgenodes(:,1) + innerstagewidth;
        
    end
    
end
function FemmProblem = slottedcommonfemmprob_radial(FemmProblem, design, ...
                            Inputs, magcornerids, Rs, coillabellocs, yokenodeids, ...
                            ymidpoint, BackIronMatInd, YokeMatInd, CoilMatInd, ...
                            YokeGroup, linktb)
% Constructs common aspects of a slotted mfemm radial flux FemmProblem
% structure
%
% Syntax
%
% FemmProblem = slottedcommonfemmprob_radial(FemmProblem, design, ...
%                             Inputs, gapedgenodes, Rs, coillabellocs, yokenodeids, ...
%                             ymidpoint, BackIronMatInd, YokeMatInd, CoilMatInd)
%
%  
% Description
%
% slottedcommonfemmprob_radial performs some final problem creation tasks
% common to all radial flux type rotary machines. These tasks include
% drawing the inner or outer air regions and linking up and adding boundary
% conditions to the simulation. This is a low-level function intended for
% use by the main radial flux machine problem creation functions rather
% than direct use.
%

    slotsperpole = design.Qs / design.Poles;

    % define the block properties of the core region
    yokeBlockProps.BlockType = FemmProblem.Materials(BackIronMatInd).Name;
    yokeBlockProps.MaxArea = Inputs.BackIronRegionMeshSize;
    yokeBlockProps.InCircuit = '';
    yokeBlockProps.InGroup = YokeGroup;

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

    % add circuits for each winding phase
    for i = 1:design.Phases
        if i == 1
            FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i), 'TotalAmps_re', Inputs.CoilCurrent);
        else
            FemmProblem = addcircuit_mfemm(FemmProblem, num2str(i));
        end
    end

    edgenodes = [];
    
    switch Inputs.StatorType
        
        case 'si'
            % single inner facing stator (magnets inside, stator outside)
            
            routerregion = [2*design.tm, 10*design.tm]; 
            
            [edgenodes(:,1), edgenodes(:,2)] = pol2cart( ...
                                                        [ 0; ...
                                                          0; ...
                                                          0; ...
                                                          0; ...
                                                          0; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2 ], ...
                                                        [ design.Rmo; ...
                                                          design.Rmo+design.g; ...
                                                          Rs + design.ty/2; ...
                                                          Rs + design.ty/2 + routerregion(1); ...
                                                          Rs + design.ty/2 + routerregion(1) + routerregion(2); ...
                                                          Rs + design.ty/2 + routerregion(1) + routerregion(2); ...
                                                          Rs + design.ty/2 + routerregion(1); ...
                                                          Rs + design.ty/2; ...
                                                          design.Rmo+design.g;
                                                          design.Rmo; ] ...
                                                        ); 

            % add the nodes to the problem
            [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, edgenodes(2:9,1), edgenodes(2:9,2));
            
            % add arcs linking the outer segments
            FemmProblem = addarcsegments_mfemm(FemmProblem, ...
                                               nodeids([2,3,4]), ...
                                               nodeids([7,6,5]), ...
                                               rad2deg(repmat(2*design.thetap,1,3)));
            
            % add segments with periodic boundaries on the outer parts
            [FemmProblem, boundind(1), boundnames{1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Stator Back Iron Periodic', 4);
            [FemmProblem, boundind(end+1), boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Stator Outer Periodic', 4);
            [FemmProblem, boundind(end+1), boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Stator Outer Periodic', 4);
            [FemmProblem, boundind(end+1), boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Air Gap Periodic', 4);

            segprops = struct('BoundaryMarker', boundnames);
            
            % bottom segs
            FemmProblem = addsegments_mfemm(FemmProblem, ...
                                            nodeids([1,2,3]), ...
                                            nodeids([2,3,4]), ...
                                            segprops(1:3));
            % top segs
            FemmProblem = addsegments_mfemm(FemmProblem, ...
                                            nodeids([8,7,6]), ...
                                            nodeids([7,6,5]), ...
                                            segprops(1:3));
            
            % bottom gap corner to bottom core corner
            FemmProblem = addsegments_mfemm(FemmProblem, nodeids(1), magcornerids(1), ...
                                               'BoundaryMarker', boundnames{4});

            % top gap corner to top core corner
            FemmProblem = addsegments_mfemm(FemmProblem, nodeids(end), magcornerids(2), ...
                                               'BoundaryMarker', boundnames{4});

            % bottom slot to edge
            FemmProblem = addarcsegments_mfemm(FemmProblem, nodeids(1), yokenodeids(1), ...
                                               rad2deg(((2*pi/design.Qs)-design.thetac)/2));

            % top slot to edge
            FemmProblem = addarcsegments_mfemm(FemmProblem, nodeids(end), yokenodeids(4), ...
                                               rad2deg(((2*pi/design.Qs)-design.thetac)/2));

            % Add block labels for the air gap
            [labelloc(1),labelloc(2)]  = pol2cart(design.thetap, design.Rmo+design.g/2);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(1).Name, ...
                                    'MaxArea', Inputs.AirGapMeshSize);

            % add a block label for the yoke and teeth
            [labelloc(1),labelloc(2)] = pol2cart(design.thetap, Rs);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                                                 yokeBlockProps);
                                        
            % Add block labels for the outer air regions
            [labelloc(1),labelloc(2)]  = pol2cart(design.thetap, Rs + design.ty/2 + routerregion(1)/2);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(1).Name, ...
                                    'MaxArea', Inputs.OuterRegionsMeshSize(1));
                                
            [labelloc(1),labelloc(2)]  = pol2cart(design.thetap, Rs + design.ty/2 + routerregion(1) + routerregion(2)/2);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(1).Name, ...
                                    'MaxArea', Inputs.OuterRegionsMeshSize(2));

                                
        case 'so'
            % single outer facing stator
            
            routerregion = [0.8*design.Ryi, 0.5*design.Ryi]; 
            
            [edgenodes(:,1), edgenodes(:,2)] = pol2cart(  ...
                                                        [ 0; ...
                                                          0; ...
                                                          0; ...
                                                          0; ...
                                                          0; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap * 2; ...
                                                          design.thetap; ...
                                                          design.thetap; ...
                                                          design.thetap ], ...
                                                        [ design.Rmi; ...
                                                          design.Rmi-design.g; ...
                                                          design.Ryi; ...
                                                          routerregion(1); ...
                                                          routerregion(2); ...
                                                          routerregion(2); ...
                                                          routerregion(1); ...
                                                          design.Ryi; ...
                                                          design.Rmi-design.g;
                                                          design.Rmi; ...
                                                          design.Ryi; ...
                                                          routerregion(1); ...
                                                          routerregion(2); ] ...
                                                        ); 

            
            
            if linktb
                
               % add the nodes to the problem
               [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, ...
                                                              edgenodes([2:5,11:13],1), ...
                                                              edgenodes([2:5,11:13],2));
                                                          
                links = [ nodeids([2,5,3,6,4,7]); ...
                          nodeids([5,2,6,3,7,4]) ];
                      
                topnodeid = nodeids(1);
                
                arcangle = design.thetap;
                
            else
                % add the nodes to the problem
                [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, ...
                                                              edgenodes(2:9,1), ...
                                                              edgenodes(2:9,2));
                                                          
                links = [ nodeids([2,3,4]); ...
                          nodeids([7,6,5]) ];
                topnodeid = nodeids(8);
                arcangle = 2*design.thetap;
            end
            
            % add arcs linking the outer segments
            FemmProblem = addarcsegments_mfemm(FemmProblem, ...
                                               links(1,:), ...
                                               links(2,:), ...
                                               rad2deg(repmat(arcangle,size(links))));
            
            if linktb
                boundnames = repmat({''}, 4,1);
                segprops = struct('BoundaryMarker', boundnames);
            else
                % add segments with periodic boundaries on the outer parts
                [FemmProblem, boundind(1), boundnames{1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Stator Back Iron Periodic', 4);
                [FemmProblem, boundind(end+1), boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Stator Outer Periodic', 4);
                [FemmProblem, boundind(end+1), boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Stator Outer Periodic', 4);
                [FemmProblem, boundind(end+1), boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'Radial Air Gap Periodic', 4);

                segprops = struct('BoundaryMarker', boundnames);

                % top segs
                FemmProblem = addsegments_mfemm(FemmProblem, ...
                                                nodeids([8,7,6]), ...
                                                nodeids([7,6,5]), ...
                                                segprops(1:3));
            end
            
            % bottom segs
            FemmProblem = addsegments_mfemm(FemmProblem, ...
                                            nodeids([1,2,3]), ...
                                            nodeids([2,3,4]), ...
                                            segprops(1:3));
                                        
            % Add block labels for the air gap
            [labelloc(1),labelloc(2)] = pol2cart(design.thetap, design.Rmi-design.g/2);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(1).Name, ...
                                    'MaxArea', Inputs.AirGapMeshSize);

            % bottom gap corner to bottom core corner
            FemmProblem = addsegments_mfemm(FemmProblem, nodeids(1), magcornerids(1), ...
                                               'BoundaryMarker', boundnames{4});

            if ~linktb
                % top gap corner to top core corner
                FemmProblem = addsegments_mfemm(FemmProblem, topnodeid, magcornerids(2), ...
                                                   'BoundaryMarker', boundnames{4});
            end
            % bottom slot to edge
            FemmProblem = addarcsegments_mfemm(FemmProblem, nodeids(1), yokenodeids(2), ...
                                               rad2deg(((2*pi/design.Qs)-design.thetac)/2));

            % top slot to edge
            FemmProblem = addarcsegments_mfemm(FemmProblem, yokenodeids(3), topnodeid, ...
                                               rad2deg(((2*pi/design.Qs)-design.thetac)/2));

            % add a block label for the yoke and teeth
            [labelloc(1),labelloc(2)] = pol2cart(design.thetap, Rs);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                                                 yokeBlockProps);

            % Add block labels for the outer air regions
            [labelloc(1),labelloc(2)]  = pol2cart(design.thetap, ...
                                            routerregion(1) + (Rs - design.ty/2 - routerregion(1))/2 );

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(1).Name, ...
                                    'MaxArea', Inputs.OuterRegionsMeshSize(1));
                                
            [labelloc(1),labelloc(2)]  = pol2cart(design.thetap, ...
                                                    routerregion(2) + (routerregion(1) - routerregion(2))/2);

            FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(labelloc(1,1), labelloc(1,2), ...
                                    'BlockType', FemmProblem.Materials(1).Name, ...
                                    'MaxArea', Inputs.OuterRegionsMeshSize(2));
                                
            if linktb
                FemmProblem.BlockLabels(end+1) = newblocklabel_mfemm(0, 0, ...
                                    'BlockType', '<No Mesh>');
            end
            
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
            error('Unrecognised StatorType option.')
                
    end

    % join up the air gaps at the top and bottom

%     % bottom right gap corner to bottom right core corner
%     FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(2), outeryokenodeids(2), ...
%                             'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);
% 
%     % top right gap corner to top right core corner
%     FemmProblem.Segments(end+1) = newsegment_mfemm(gapcornernodeids(3), outeryokenodeids(3), ...
%                             'BoundaryMarker', FemmProblem.BoundaryProps(boundind(3)).Name);
% 
%     % close the tooth and air boundaries to complete the core
%     FemmProblem = addsegments_mfemm(FemmProblem, yokenodeids(i,1:4), outeryokenodeids(1:4));  

%         corexpos = corexpos + innerstagewidth;


    % add block labels for the coils
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

            FemmProblem = addblocklabel_mfemm( FemmProblem, ...
                                               coillabellocs(row,1), ...
                                               coillabellocs(row,2), ...
                                               coilBlockProps);

%                 row = row + 1;
            if row+(Inputs.NSlots*Inputs.NWindingLayers) <= size(coillabellocs,1)
                FemmProblem = addblocklabel_mfemm( FemmProblem, ...
                                                   coillabellocs(row+(Inputs.NSlots*Inputs.NWindingLayers),1), ...
                                                   coillabellocs(row,2), ...
                                                   coilBlockProps );
            end
            
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
        
    
end
function [FemmProblem, outernodes, coillabellocs] = radialfluxstatorhalf2dfemmprob(...
    slots, Poles, thetapole, thetacoil, thetashoegap, ryoke, rcoil, rshoebase, rshoegap, roffset, side, varargin)
% draw internal parts of half a slotted stator
%
% Syntax
%
% [FemmProblem, outernodes, coillabellocs] = ...
%       radialfluxstatorhalf2dfemmprob(slots, Poles, thetapole, thetacoil, thetashoegap, ...
%       ryoke, rcoil, rshoebase, rshoegap, coillayers, side, varargin)
%
%
% Input
%

    Inputs.NWindingLayers = 1;
    Inputs.FemmProblem = newproblem_mfemm('planar');
    Inputs.ToothMaterial = 1;
    Inputs.ToothRegionMeshSize = -1;
    Inputs.ShoeGapMaterial = 1;
    Inputs.ShoeGapRegionMeshSize = -1;
    Inputs.ShoeGroup = 1;
    Inputs.SlotPositions = [];
    Inputs.NSlots = [];
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % copy over the FemmProblem
    FemmProblem = Inputs.FemmProblem;
    
    % remove the odl one to save memory
    Inputs = rmfield(Inputs, 'FemmProblem');
    
    % check an integer number of machine Poles and slots (can't have half a
    % slot or a pole in a machine
    slots = round(slots);
    Poles = round(Poles);

    slotsperpole = slots / Poles;

    % The user can either supply a set of slot positions (useful for
    % inductance simulations), or these will be calculated to fill two
    % Poles of the machine
    if isempty(Inputs.SlotPositions)
        
        if isempty(Inputs.NSlots)
            
            toothpos = 0:1/slotsperpole:2;

        else
            
            % draw a set number of slots
            toothpos = linspace(0, Inputs.NSlots*1/slotsperpole, Inputs.NSlots+1);

        end
        
        % distribute the slots between the denormalised tooth positions, and
        % add a full slot separation as the first slot is drawn split along the
        % line y = 0
        slotpos = toothpos(1:end-1)*thetapole + (thetapole / slotsperpole / 2);

    else
        slotpos = Inputs.SlotPositions;
    end
    
    elcount = elementcount_mfemm(FemmProblem);
    
    % make a single slot
    [nodes, links, cornernodes, shoelabelloc, shoegaplabelloc, firstcoillabellocs] = ...
            internalslotnodelinks(thetacoil, thetashoegap, ryoke/2, rcoil, rshoebase, rshoegap, Inputs.NWindingLayers, Inputs.Tol);

    coillabellocs = [];
    
    % we flip the node positions if we are drawing an internally facing
    % stator
    if strcmp(side, 'i')
        
        nodes(:,1) = -nodes(:,1);
        
        if ~isempty(shoelabelloc)
            shoelabelloc(:,1) = -shoelabelloc(:,1);
        end
        
        if ~isempty(shoegaplabelloc)
            shoegaplabelloc(:,1) = -shoegaplabelloc(:,1);
        end
        
        if ~isempty(firstcoillabellocs)
            firstcoillabellocs(:,1) = -firstcoillabellocs(:,1);
        end
        
        % rearrange the corner nodes to preserve clockwise ordering
        % starting from bottom left
        cornernodes = [cornernodes(2), cornernodes(1), cornernodes(4), cornernodes(3)]; 
    end
    
    % add the specified offset in the radial direction
    nodes(:,1) = nodes(:,1) + roffset;

    if ~isempty(shoelabelloc)
        shoelabelloc(:,1) = shoelabelloc(:,1) + roffset;
    end

    if ~isempty(shoegaplabelloc)
        shoegaplabelloc(:,1) = shoegaplabelloc(:,1) + roffset;
    end

    if ~isempty(firstcoillabellocs)
        firstcoillabellocs(:,1) = firstcoillabellocs(:,1) + roffset;
    end
    
	% get the vertical links by finding those links where the difference in
    % y coordinates of the link nodes is not zero, these links must be made
    % into arc segments
    isvertlinks = abs(diff( [nodes(links(:,1)+1,1), nodes(links(:,2)+1,1)], 1, 2 )) < Inputs.Tol;
    vertlinks = links(isvertlinks,:);
    angles = diff( [nodes(vertlinks(:,1)+1,2), nodes(vertlinks(:,2)+1,2)], 1, 2);
    
    % get the horizontal links, these will be segments
    horizlinks = links(~isvertlinks,:);
    
    % correct vertical links which are in the wrong direction
    for i = 1:size(vertlinks,1)
        if angles(i) < 0
            vertlinks(i,:) = fliplr(vertlinks(i,:));
            angles(i) = abs(angles(i));
        end
    end
    % convert angles to degrees
    angles = rad2deg(angles);
    vertlinks = fliplr(vertlinks);
    
    % transform the node locations to convert the rectangulr region to the
    % desired arced region 
    [nodes(:,1), nodes(:,2)] = pol2cart(nodes(:,2), nodes(:,1));
    if ~isempty(shoelabelloc)
        [shoelabelloc(:,1), shoelabelloc(:,2)] = pol2cart(shoelabelloc(:,2),shoelabelloc(:,1));
    end
    if ~isempty(shoegaplabelloc)
        [shoegaplabelloc(:,1), shoegaplabelloc(:,2)] = pol2cart(shoegaplabelloc(:,2),shoegaplabelloc(:,1));
    end
    [firstcoillabellocs(:,1), firstcoillabellocs(:,2)] = pol2cart(firstcoillabellocs(:,2),firstcoillabellocs(:,1));

    % store the nodes at the bottom of all the slots
    bottomnodes = cornernodes([4,3]) + elcount.NNodes;  
    lastslotcornernodes = cornernodes + elcount.NNodes;    
    
    % draw the first slot linking it to the bottom of the domain
    thisslotsnodes = nodes;
        
    rotM = [cos(slotpos(1))  sin(slotpos(1)); 
            sin(slotpos(1)) -cos(slotpos(1))];
            
    % move in the y direction to the first slot position
    thisslotsnodes = thisslotsnodes * rotM;

    if ~isempty(firstcoillabellocs)
        
        thiscoillabellocs = firstcoillabellocs * rotM;
        
        % get the coil label location
        coillabellocs = [ coillabellocs; ...
                          thiscoillabellocs ];
                      
    end
        
    thisslotvertlinks = vertlinks + elcount.NNodes;
    thisslothorizlinks = horizlinks + elcount.NNodes;
    
    [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, thisslotsnodes(:,1), thisslotsnodes(:,2));
        
    [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, ...
                                               thisslothorizlinks(:,1), ...
                                               thisslothorizlinks(:,2));
                                           
    [FemmProblem, seginds] = addarcsegments_mfemm(FemmProblem, ...
                                                  thisslotvertlinks(:,1), ...
                                                  thisslotvertlinks(:,2), ...
                                                  angles, ...
                                                  'MaxSegDegrees', angles ./ 20);
        
    if ~isempty(shoelabelloc) 
            
        thisshoelabelloc = shoelabelloc * rotM;
        
        % add the block labels for the shoes and gap
        FemmProblem = addblocklabel_mfemm(FemmProblem, thisshoelabelloc(1,1), thisshoelabelloc(1,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.ToothMaterial).Name, ...
                                            'MaxArea', Inputs.ToothRegionMeshSize, ...
                                            'InGroup', Inputs.ShoeGroup);

        FemmProblem = addblocklabel_mfemm(FemmProblem, thisshoelabelloc(2,1), thisshoelabelloc(2,2), ...
                                            'BlockType', FemmProblem.Materials(Inputs.ToothMaterial).Name, ...
                                            'MaxArea', Inputs.ToothRegionMeshSize, ...
                                            'InGroup', Inputs.ShoeGroup);
    
    end
    
    if ~isempty(shoegaplabelloc)
        
        thisshoegaplabelloc = shoegaplabelloc * rotM;
        
        for ind = 1:size(thisshoegaplabelloc,1)
            FemmProblem = addblocklabel_mfemm(FemmProblem, thisshoegaplabelloc(ind,1), thisshoegaplabelloc(ind,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
                                                'MaxArea', Inputs.ShoeGapRegionMeshSize);
        end
        
    end
    
    
    % draw the rest of the slots, making copies of the original slot's
    % nodes and links, adding in the labels, and linking them together
    for i = 2:numel(slotpos)
        
        thisslotsnodes = nodes;
        
        rotM = [cos(slotpos(i))  sin(slotpos(i)); 
                sin(slotpos(i)) -cos(slotpos(i))]; 
        
        thisslotsnodes = thisslotsnodes * rotM;
        
        thisslotvertlinks = vertlinks + numel(FemmProblem.Nodes);
        thisslothorizlinks = horizlinks + numel(FemmProblem.Nodes);
        
        thisslotcornernodes = cornernodes + numel(FemmProblem.Nodes);
        
        [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, ...
                                                          thisslotsnodes(:,1), ...
                                                          thisslotsnodes(:,2));
        
        [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, ...
                                                   thisslothorizlinks(:,1), ...
                                                   thisslothorizlinks(:,2));
        
        [FemmProblem, seginds] = addarcsegments_mfemm(FemmProblem, ...
                                                      thisslotvertlinks(:,1), ...
                                                      thisslotvertlinks(:,2), ...
                                                      angles, ...
                                                      'MaxSegDegrees', angles ./ 20 );
                                               
        % add the block labels for the shoes and gap
        if ~isempty(shoelabelloc)
            
            thisshoelabelloc = shoelabelloc * rotM;
        
            % add the block labels for the shoes and gap
            FemmProblem = addblocklabel_mfemm(FemmProblem, thisshoelabelloc(1,1), thisshoelabelloc(1,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.ToothMaterial).Name, ...
                                                'MaxArea', Inputs.ToothRegionMeshSize, ...
                                                'InGroup', Inputs.ShoeGroup);

            FemmProblem = addblocklabel_mfemm(FemmProblem, thisshoelabelloc(2,1), thisshoelabelloc(2,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.ToothMaterial).Name, ...
                                                'MaxArea', Inputs.ToothRegionMeshSize, ...
                                                'InGroup', Inputs.ShoeGroup);
    
        end
        
        if ~isempty(shoegaplabelloc)
            
            thisshoegaplabelloc = shoegaplabelloc * rotM;
            
            for ind = 1:size(thisshoegaplabelloc,1)
                
                FemmProblem = addblocklabel_mfemm(FemmProblem, thisshoegaplabelloc(ind,1), thisshoegaplabelloc(ind,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
                                                'MaxArea', Inputs.ShoeGapRegionMeshSize);

            end
                                        
        end
        
        % now link the slot to the slot below
        if side == 'i'
            
            slotangles = chordpoints2angle(FemmProblem.Nodes(thisslotcornernodes(4)+1).Coords(1), ...
                                       FemmProblem.Nodes(thisslotcornernodes(4)+1).Coords(2), ...
                                       FemmProblem.Nodes(lastslotcornernodes(1)+1).Coords(1), ...
                                       FemmProblem.Nodes(lastslotcornernodes(1)+1).Coords(2) );
                                   
            [FemmProblem, seginds] = ...
                addarcsegments_mfemm( FemmProblem, ...
                                      lastslotcornernodes(1), ...
                                      thisslotcornernodes(4), ...
                                      slotangles, ...
                                      'MaxSegDegrees', slotangles ./ 20 );
                                  
        else
%             [FemmProblem, seginds] = ...
%                 addsegments_mfemm(FemmProblem, thisslotcornernodes(2), lastslottopnodes(1));
%             
            slotangles = chordpoints2angle(FemmProblem.Nodes(thisslotcornernodes(3)+1).Coords(1), ...
                                       FemmProblem.Nodes(thisslotcornernodes(3)+1).Coords(2), ...
                                       FemmProblem.Nodes(lastslotcornernodes(2)+1).Coords(1), ...
                                       FemmProblem.Nodes(lastslotcornernodes(2)+1).Coords(2) );
                                   
            [FemmProblem, seginds] = ...
                addarcsegments_mfemm( FemmProblem, ...
                                      lastslotcornernodes(2), ...
                                      thisslotcornernodes(3), ...
                                      slotangles, ...
                                      'MaxSegDegrees', slotangles ./ 20 );

        end
        
        % store the top nodes of the slot for the next loop
        lastslotcornernodes = thisslotcornernodes;     
        
        thiscoillabellocs = firstcoillabellocs * rotM;
        
        % get the coil label location
        coillabellocs = [ coillabellocs; ...
                          thiscoillabellocs ];
                      
%                           coillabellocs(1:Inputs.NWindingLayers,1), ...
%                           coillabellocs(1:Inputs.NWindingLayers,2) + slotpos(i) - slotpos(1)];
        
    end
    
    % we will return the outer corner node ids for later use
    outernodes = [bottomnodes, lastslotcornernodes([2,1])];
    
end

function angle = chordpoints2angle(x1,y1,x2,y2)

     chordlen = hypot(x2-x1, y2-y1);
     
     R = magn([x1,y1]);
     
     angle = chordlen/R;
     
     angle = rad2deg(angle);

end


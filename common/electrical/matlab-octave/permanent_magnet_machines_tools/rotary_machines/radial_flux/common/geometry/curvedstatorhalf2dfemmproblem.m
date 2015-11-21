function [FemmProblem, outernodes, coillabellocs, slotinfo] = curvedstatorhalf2dfemmproblem(...
    nslots, thetaslot, thetacoil, thetashoegap, ryoke, rcoil, rshoebase, rshoegap, roffset, side, varargin)
% draw internal parts of a slotted stator for a radial flux machine
%
% Syntax
%
% [FemmProblem, outernodes, coillabellocs] = ...
%       curvedstatorhalf2dfemmproblem( nslots, thetaslot, thetacoil, ...
%                                       thetashoegap, ryoke, rcoil, rshoebase, ...
%                                       rshoegap, coillayers, side )
%
%
% [FemmProblem, outernodes, coillabellocs] = ...
%       curvedstatorhalf2dfemmproblem( ..., 'Parameter', Value )
%
% Description
%
% curvedstatorhalf2dfemmproblem creates or adds to an mfemm FemmProblem
% structure, the internal parts of a radial flux slotted stator. The
% 'internal parts' refer to the segments, nodes and labels making up the
% teeth and coils, but not the outer part of the yoke of the machine, i.e.
% the stator drawing is not closed at its edges. It is intended for
% curvedstatorhalf2dfemmproblem to be used by higher level drawing
% functions which are drawing the full radial flux machine geometry. 
% 
%
% Input
%
%  nslots - total number of slots to be drawn
%
%  thetaslot - angular displacement between slots
%
%  thetacoil - 
%
%  thetashoegap - angular span of slot opening
%
%  ryoke - radial thickness of the sttor yoke
%
%  rcoil - radial height of the filled slot
%
%  rshoebase - 
%
%  rshoegap - 
%
%  roffset - 
%
%  In addition, a number of other optional arguments can be supplied as
%  parameter-value pairs. These options and their behaviour are as follows:
%
%  'FemmProblem'
%  'NWindingLayers'
%  'SlotPositions'
%  'NSlots'
%  'Tol'
%  'ShoeGapMaterial'
%  'ShoeGapRegionMeshSize'
%  'CoilBaseFraction'
%  'ShoeCurveControlFrac'
%  'SplitX'
%  'CoilInsulationThickness'
%  'DrawCoilInsulation'
%
% Output
%
% 

    Inputs.NWindingLayers = 1;
    Inputs.FemmProblem = newproblem_mfemm('planar');
    Inputs.ShoeGapMaterial = 1;
    Inputs.ShoeGapRegionMeshSize = -1;
%     Inputs.ArmatureIronGroup = 1;
    Inputs.Tol = 1e-5;
    Inputs.CoilBaseFraction = 0.05;
    Inputs.ShoeCurveControlFrac = 0.5;
    Inputs.SplitX = false;
    Inputs.DrawCoilInsulation = false;
    Inputs.CoilInsulationThickness = 0;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % copy over the FemmProblem
    FemmProblem = Inputs.FemmProblem;
    
    % remove the old one to save memory
    Inputs = rmfield(Inputs, 'FemmProblem');

    slotpos = linspace(0, (nslots-1)*thetaslot, nslots) + (thetaslot / 2);
    
    elcount = elementcount_mfemm(FemmProblem);
    
    if Inputs.DrawCoilInsulation
        insulationthickness = Inputs.CoilInsulationThickness;
    else
        insulationthickness = 0;
    end
    
    % make a single slot
    [nodes, links, slotinfo] = internalslotnodelinks( thetacoil, thetashoegap, ryoke/2, rcoil, ...
                                   rshoebase, rshoegap, Inputs.NWindingLayers, Inputs.Tol, ...
                                   'CoilBaseFraction', Inputs.CoilBaseFraction, ...
                                   'InsulationThickness', insulationthickness, ...
                                   'ShoeCurveControlFrac', Inputs.ShoeCurveControlFrac, ...
                                   'YScale', roffset + ryoke + rcoil/2  );

    links = [ links, ismember(1:size(links,1),slotinfo.vertlinkinds)', ismember(1:size(links,1),slotinfo.toothlinkinds)' ];
    
    coillabellocs = [];
    
    % we flip the node positions if we are drawing an internally facing
    % stator
    if strcmp(side, 'i')
        
        nodes(:,1) = -nodes(:,1);
        
        if ~isempty(slotinfo.shoegaplabelloc)
            slotinfo.shoegaplabelloc(:,1) = -slotinfo.shoegaplabelloc(:,1);
        end
        
        if ~isempty(slotinfo.coillabelloc)
            slotinfo.coillabelloc(:,1) = -slotinfo.coillabelloc(:,1);
        end
        
        if ~isempty (slotinfo.inslabelloc)
            slotinfo.inslabelloc(:,1) = -slotinfo.inslabelloc(:,1);
        end
        
        % rearrange the corner nodes to preserve clockwise ordering
        % starting from bottom left
        slotinfo.cornernodes = [slotinfo.cornernodes(2), slotinfo.cornernodes(1), slotinfo.cornernodes(4), slotinfo.cornernodes(3)]; 
    end
    
    % add the specified offset in the radial direction
    nodes(:,1) = nodes(:,1) + roffset;

    if ~isempty(slotinfo.shoegaplabelloc)
        slotinfo.shoegaplabelloc(:,1) = slotinfo.shoegaplabelloc(:,1) + roffset;
    end

    if ~isempty(slotinfo.coillabelloc)
        slotinfo.coillabelloc(:,1) = slotinfo.coillabelloc(:,1) + roffset;
    end
    
    if ~isempty(slotinfo.inslabelloc)
        slotinfo.inslabelloc(:,1) = slotinfo.inslabelloc(:,1) + roffset;
    end
    
	% convert vertical links to arc segments
    vertlinks = links(slotinfo.vertlinkinds,1:2);
    angles = diff( [nodes(vertlinks(:,1)+1,2), nodes(vertlinks(:,2)+1,2)], 1, 2);
    
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
    links(slotinfo.vertlinkinds, 1:2) = vertlinks;
    
    % transform the node locations to convert the rectangular region to the
    % desired arced region 
    [nodes(:,1), nodes(:,2)] = pol2cart(nodes(:,2), nodes(:,1));
    if ~isempty(slotinfo.shoegaplabelloc)
        [slotinfo.shoegaplabelloc(:,1), slotinfo.shoegaplabelloc(:,2)] = pol2cart(slotinfo.shoegaplabelloc(:,2),slotinfo.shoegaplabelloc(:,1));
    end
    
    [slotinfo.coillabelloc(:,1), slotinfo.coillabelloc(:,2)] = pol2cart(slotinfo.coillabelloc(:,2),slotinfo.coillabelloc(:,1));
    
    if ~isempty(slotinfo.inslabelloc)
        [slotinfo.inslabelloc(:,1), slotinfo.inslabelloc(:,2)] = pol2cart(slotinfo.inslabelloc(:,2),slotinfo.inslabelloc(:,1));
    end
    
    % store the nodes at the bottom of all the slots
    bottomnodes = slotinfo.cornernodes([4,3]) + elcount.NNodes;  
    lastslotcornernodes = slotinfo.cornernodes + elcount.NNodes;    
    
    % draw the first slot linking it to the bottom of the domain
    thisslotsnodes = nodes;
        
    rotM = [cos(slotpos(1))  sin(slotpos(1)); 
            sin(slotpos(1)) -cos(slotpos(1))];
            
    % move in the y direction to the first slot position
    thisslotsnodes = thisslotsnodes * rotM;
    
    originslabellocs = slotinfo.inslabelloc;
    
    if ~isempty (originslabellocs)
        slotinfo.inslabelloc = originslabellocs * rotM;
    end
    
    if ~isempty(slotinfo.coillabelloc)
        
        thiscoillabellocs = slotinfo.coillabelloc * rotM;
        
        % get the coil label location
        coillabellocs = [ coillabellocs; ...
                          thiscoillabellocs ];

    end
    
    % add a new group number for the stator iron
    if ~isfield (FemmProblem.Groups, 'StatorIronOutline')
        FemmProblem = addgroup_mfemm (FemmProblem, 'StatorIronOutline');
    end
    statorirongp = getgroupnumber_mfemm (FemmProblem, 'StatorIronOutline');

    % TODO: put iron and insulation in different groups
%    ironinds = setdiff(1:size(links,1),inslinkinds);
    thisslotlinks = [links(:,1:2) + numel(FemmProblem.Nodes), links(:,3:end)];
    
    [FemmProblem, ~, ~] = addnodes_mfemm (FemmProblem, thisslotsnodes(:,1), thisslotsnodes(:,2));
    
    vcount = 1;
    for ii = 1:size(thisslotlinks,1)
        
        if thisslotlinks (ii,3)
            % vertical link that must be converted to an arc
            if thisslotlinks (ii,4)
                
               [FemmProblem, ~] = addarcsegments_mfemm (FemmProblem, ...
                                                  thisslotlinks(ii,1), ...
                                                  thisslotlinks(ii,2), ...
                                                  angles(vcount), ...
                                                  'MaxSegDegrees', angles(vcount) ./ max (4, ceil (angles(vcount)/1)), ...
                                                  'InGroup', statorirongp);
            else
               [FemmProblem, ~] = addarcsegments_mfemm (FemmProblem, ...
                                                  thisslotlinks(ii,1), ...
                                                  thisslotlinks(ii,2), ...
                                                  angles(vcount), ...
                                                  'MaxSegDegrees', angles(vcount) ./ max (4, ceil (angles(vcount)/1))); 
            end
            vcount = vcount + 1;
        else
            % horizontal link
            if thisslotlinks (ii,4)
                [FemmProblem, ~] = addsegments_mfemm (FemmProblem, ...
                                                      thisslotlinks(ii,1), ...
                                                      thisslotlinks(ii,2), ...
                                                      'InGroup', statorirongp);
            else
                [FemmProblem, ~] = addsegments_mfemm (FemmProblem, ...
                                                      thisslotlinks(ii,1), ...
                                                      thisslotlinks(ii,2));
            end
        end
        
    end
          
    if ~isempty (slotinfo.shoegaplabelloc)
        
        thisshoegaplabelloc = slotinfo.shoegaplabelloc * rotM;
        
        for ind = 1:size(thisshoegaplabelloc,1)
            FemmProblem = addblocklabel_mfemm (FemmProblem, thisshoegaplabelloc(ind,1), thisshoegaplabelloc(ind,2), ...
                                                'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
                                                'MaxArea', Inputs.ShoeGapRegionMeshSize );
        end
        
    end
    
    % draw the rest of the slots, making copies of the original slot's
    % nodes and links, adding in the labels, and linking them together
    for i = 2:numel(slotpos)
    
        slotinsgp = getgroupnumber_mfemm (FemmProblem, ['StatorSlot_', int2str(i)]);
        
        thisslotsnodes = nodes;
        
        rotM = [cos(slotpos(i))  sin(slotpos(i)); 
                sin(slotpos(i)) -cos(slotpos(i))]; 
        
        thisslotsnodes = thisslotsnodes * rotM;
        
        if ~isempty (slotinfo.inslabelloc)
            slotinfo.inslabelloc = [ slotinfo.inslabelloc; originslabellocs * rotM ];
        end
        
        thisslotlinks = [links(:,1:2) + numel(FemmProblem.Nodes), links(:,3:end)];
        
        thisslotcornernodes = slotinfo.cornernodes + numel(FemmProblem.Nodes);
        
        [FemmProblem] = addnodes_mfemm (FemmProblem, ...
                                        thisslotsnodes(:,1), ...
                                        thisslotsnodes(:,2));
        
        vcount = 1;
        for ii = 1:size(thisslotlinks,1)

            if thisslotlinks (ii,3)
                % vertical link that must be converted to an arc
                if thisslotlinks (ii,4)

                   [FemmProblem, ~] = addarcsegments_mfemm(FemmProblem, ...
                                                      thisslotlinks(ii,1), ...
                                                      thisslotlinks(ii,2), ...
                                                      angles(vcount), ...
                                                      'MaxSegDegrees', angles(vcount) ./ max (4, ceil (angles(vcount)/1)), ...
                                                      'InGroup', statorirongp);
                else
                   [FemmProblem, ~] = addarcsegments_mfemm(FemmProblem, ...
                                                      thisslotlinks(ii,1), ...
                                                      thisslotlinks(ii,2), ...
                                                      angles(vcount), ...
                                                      'MaxSegDegrees', angles(vcount) ./ max (4, ceil (angles(vcount)/1))); 
                end
                vcount = vcount + 1;
            else
                % horizontal link
                if thisslotlinks (ii,4)
                    [FemmProblem, ~] = addsegments_mfemm (FemmProblem, ...
                                                          thisslotlinks(ii,1), ...
                                                          thisslotlinks(ii,2), ...
                                                          'InGroup', statorirongp);
                else
                    [FemmProblem, ~] = addsegments_mfemm (FemmProblem, ...
                                                          thisslotlinks(ii,1), ...
                                                          thisslotlinks(ii,2));
                end
            end

        end                              
        
        if ~isempty(slotinfo.shoegaplabelloc)
            
            thisshoegaplabelloc = slotinfo.shoegaplabelloc * rotM;
            
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
                                      'MaxSegDegrees', slotangles ./ 20, ...
                                      'InGroup', statorirongp );
                                  
        else
            
            slotangles = chordpoints2angle(FemmProblem.Nodes(thisslotcornernodes(3)+1).Coords(1), ...
                                       FemmProblem.Nodes(thisslotcornernodes(3)+1).Coords(2), ...
                                       FemmProblem.Nodes(lastslotcornernodes(2)+1).Coords(1), ...
                                       FemmProblem.Nodes(lastslotcornernodes(2)+1).Coords(2) );
                                   
            [FemmProblem, seginds] = ...
                addarcsegments_mfemm( FemmProblem, ...
                                      lastslotcornernodes(2), ...
                                      thisslotcornernodes(3), ...
                                      slotangles, ...
                                      'MaxSegDegrees', slotangles ./ 20, ...
                                      'InGroup', statorirongp );

        end
        
        % store the top nodes of the slot for the next loop
        lastslotcornernodes = thisslotcornernodes;     
        
        thiscoillabellocs = slotinfo.coillabelloc * rotM;
        
        % get the coil label location
        coillabellocs = [ coillabellocs; ...
                          thiscoillabellocs ];
                           
    end
    
    % we will return the outer corner node ids for later use
    outernodes = [bottomnodes, lastslotcornernodes([2,1])];
    
    
    
    % ---------------    Nested functions   -------------------- %
    
    
%    % nested function for drawing single slots
%    function drawslot ()
%    
%        vcount = 1;
%        for i = 1:size(thisslotlinks,1)
%            
%            if thisslotlinks (i,3)
%                % vertical link that must be converted to an arc
%                if thisslotlinks (i,4)
%                    
%                   [FemmProblem, ~] = addarcsegments_mfemm (FemmProblem, ...
%                                                      thisslotlinks(i,1), ...
%                                                      thisslotlinks(i,2), ...
%                                                      angles(vcount), ...
%                                                      'MaxSegDegrees', angles(vcount) ./ 20, ...
%                                                      'InGroup', statorirongp);
%                else
%                   [FemmProblem, ~] = addarcsegments_mfemm (FemmProblem, ...
%                                                      thisslotlinks(i,1), ...
%                                                      thisslotlinks(i,2), ...
%                                                      angles(vcount), ...
%                                                      'MaxSegDegrees', angles(vcount) ./ 20); 
%                end
%                vcount = vcount + 1;
%            else
%                % horizontal link
%                if thisslotlinks (i,4)
%                    [FemmProblem, ~] = addsegments_mfemm (FemmProblem, ...
%                                                          thisslotlinks(i,1), ...
%                                                          thisslotlinks(i,2), ...
%                                                          'InGroup', statorirongp);
%                else
%                    [FemmProblem, ~] = addsegments_mfemm (FemmProblem, ...
%                                                          thisslotlinks(i,1), ...
%                                                          thisslotlinks(i,2));
%                end
%            end
%            
%        end
%              
%        if ~isempty (slotinfo.shoegaplabelloc)
%            
%            thisshoegaplabelloc = slotinfp.shoegaplabelloc * rotM;
%            
%            for ind = 1:size(thisshoegaplabelloc,1)
%                FemmProblem = addblocklabel_mfemm (FemmProblem, thisshoegaplabelloc(ind,1), thisshoegaplabelloc(ind,2), ...
%                                                    'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
%                                                    'MaxArea', Inputs.ShoeGapRegionMeshSize );
%            end
%            
%        end      
%    
%    end
    
end


function angle = chordpoints2angle(x1,y1,x2,y2)

     chordlen = hypot(x2-x1, y2-y1);
     
     R = magn([x1,y1]);
     
     angle = chordlen/R;
     
     angle = rad2deg(angle);

end


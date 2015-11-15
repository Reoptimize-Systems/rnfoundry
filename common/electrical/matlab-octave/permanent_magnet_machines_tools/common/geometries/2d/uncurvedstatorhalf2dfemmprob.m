function [FemmProblem, outernodes, coillabellocs, slotinfo] = uncurvedstatorhalf2dfemmprob (nslots, yslot, ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, xoffset, side, varargin)
% draw internal parts of half a slotted axial flux stator
%
% Syntax
%
% [FemmProblem, outernodes, slotinfo.coillabellocs] = ...
%       axialfluxstatorhalf2dfemmprob(slots, Poles, ypole, ycoil, yshoegap, ...
%       xyoke, xcoil, xshoebase, xshoegap, xoffset, side, varargin)
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
    Inputs.Tol = 1e-5;
    Inputs.CoilBaseFraction = 0.05;
    Inputs.ShoeCurveControlFrac = 0.5;
    Inputs.SplitX = false;
    Inputs.DrawCoilInsulation = false;
    Inputs.CoilInsulationThickness = 0;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % copy over the FemmProblem
    FemmProblem = Inputs.FemmProblem;
    
    % remove the odl one to save memory
    Inputs = rmfield(Inputs, 'FemmProblem');
    
    % check an integer number of machine slots (can't have half a slot or a
    % pole in a machine
    if ~isint2eps (nslots)
        error ('Number of slots must be an integer.')
    end
    
    if Inputs.DrawCoilInsulation
        insulationthickness = Inputs.CoilInsulationThickness;
    else
        insulationthickness = 0;
    end

    % Calculate slot positions. We shift the slots up by half a slot width
    % as the internalslotnodelinks function draws the slot symmetric about
    % the line y = 0
    slotpos = linspace(0, (nslots-1)*yslot, nslots) + (yslot / 2);
    
    % make a single slot
    [nodes, links, slotinfo] = ...
            internalslotnodelinks(ycoil, yshoegap, xyoke/2, xcoil, xshoebase, xshoegap, ...
                                  Inputs.NWindingLayers, ...
                                  Inputs.Tol, ...
                                  'CoilBaseFraction', Inputs.CoilBaseFraction, ...
                                  'InsulationThickness', insulationthickness, ...
                                  'ShoeCurveControlFrac', Inputs.ShoeCurveControlFrac );
                        
	coillabellocs = [];
    
    % we flip the node positions if we are drawing the left hand side
    if side == 'l'
        
        nodes(:,1) = -nodes(:,1);

        if ~isempty(slotinfo.shoegaplabelloc)
            slotinfo.shoegaplabelloc(:,1) = -slotinfo.shoegaplabelloc(:,1);
        end
        
        if ~isempty(slotinfo.coillabelloc)
            slotinfo.coillabelloc(:,1) = -slotinfo.coillabelloc(:,1);
        end
        
        % rearrange the corner nodes to preserve clockwise ordering
        % starting from bottom left
        slotinfo.cornernodes = [slotinfo.cornernodes(2), slotinfo.cornernodes(1), slotinfo.cornernodes(4), slotinfo.cornernodes(3)];
        
    end
    
    % add the specified offset in the x direction
    nodes(:,1) = nodes(:,1) + xoffset;

    if ~isempty(slotinfo.shoegaplabelloc)
        slotinfo.shoegaplabelloc(:,1) = slotinfo.shoegaplabelloc(:,1) + xoffset;
    end

    if ~isempty(slotinfo.coillabelloc)
        slotinfo.coillabelloc(:,1) = slotinfo.coillabelloc(:,1) + xoffset;
    end
        
    elcount = elementcount_mfemm(FemmProblem);
    
    % store the nodes at the bottom of all the slots
    bottomnodes = slotinfo.cornernodes(1:2) + elcount.NNodes;  
    lastslottopnodes = slotinfo.cornernodes(3:4) + elcount.NNodes;    
    
    % draw the first slot linking it to the bottom of the domain
    thisslotsnodes = nodes;
        
    % move in the y direction to the first slot position
    thisslotsnodes(:,2) = thisslotsnodes(:,2) + slotpos(1);

    if ~isempty(slotinfo.coillabelloc)
        slotinfo.coillabelloc(:,2) = slotinfo.coillabelloc(:,2) + slotpos(1);
    end
        
    thisslotlinks = links + elcount.NNodes;
    
    [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, thisslotsnodes(:,1), thisslotsnodes(:,2));
        
    [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, thisslotlinks(:,1), thisslotlinks(:,2));
    
    if ~isempty(slotinfo.shoegaplabelloc)
        FemmProblem = addblocklabel_mfemm(FemmProblem, slotinfo.shoegaplabelloc(1,1), slotinfo.shoegaplabelloc(1,2) + slotpos(1), ...
                                            'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
                                            'MaxArea', Inputs.ShoeGapRegionMeshSize);
    end
    
    % get the coil label location
    coillabellocs = [ slotinfo.coillabelloc(1:Inputs.NWindingLayers,1), ...
                      slotinfo.coillabelloc(1:Inputs.NWindingLayers,2) + slotpos(1) ];
    
    % draw the rest of the slots, making copies of the original slot's
    % nodes and links, adding in the labels, and linking them together
    for i = 2:numel(slotpos)
        
        thisslotsnodes = nodes;
        
        thisslotsnodes(:,2) = thisslotsnodes(:,2) + slotpos(i);
        
        thisslotlinks = links + numel(FemmProblem.Nodes);
        
        thisslotcornernodes = slotinfo.cornernodes + numel(FemmProblem.Nodes);
        
        [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, thisslotsnodes(:,1), thisslotsnodes(:,2));
        
        [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, thisslotlinks(:,1), thisslotlinks(:,2));
        
        if ~isempty(slotinfo.shoegaplabelloc)
            
            FemmProblem = addblocklabel_mfemm(FemmProblem, slotinfo.shoegaplabelloc(1,1), slotinfo.shoegaplabelloc(1,2) + slotpos(i), ...
                                            'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
                                            'MaxArea', Inputs.ShoeGapRegionMeshSize);
                                        
        end
        
        % now link the slot to the slot below
        if side == 'l'
            [FemmProblem, seginds] = ...
                addsegments_mfemm(FemmProblem, thisslotcornernodes(1), lastslottopnodes(2));            
        else
            [FemmProblem, seginds] = ...
                addsegments_mfemm(FemmProblem, thisslotcornernodes(2), lastslottopnodes(1));
        end
        
        % store the top nodes of the slot for the next loop
        lastslottopnodes = thisslotcornernodes(3:4);     
        
        % get the coil label location
        coillabellocs = [ coillabellocs; ...
                          slotinfo.coillabelloc(1:Inputs.NWindingLayers,1), ...
                          slotinfo.coillabelloc(1:Inputs.NWindingLayers,2) + slotpos(i) ];
        
    end
    
    % we will return the outer corner node ids for later use
    outernodes = [bottomnodes, lastslottopnodes];
    
end


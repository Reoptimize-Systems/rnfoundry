function [FemmProblem, outernodes, coillabellocs] = axialfluxstatorhalf2dfemmprob(slots, Poles, ypole, ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, xoffset, side, varargin)
% draw internal parts of half a slotted stator
%
% Syntax
%
% [FemmProblem, outernodes, coillabellocs] = ...
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
        slotpos = toothpos(1:end-1)*ypole + (ypole / slotsperpole / 2);

    else
        slotpos = Inputs.SlotPositions;
    end
    
    % make a single slot
    [nodes, links, cornernodes, shoegaplabelloc, coillabellocs] = ...
            internalslotnodelinks(ycoil, yshoegap, xyoke/2, xcoil, xshoebase, xshoegap, Inputs.NWindingLayers, Inputs.Tol);
    
    % we flip the node positions if we are drawing the left hand side
    if side == 'l'
        
        nodes(:,1) = -nodes(:,1);

        if ~isempty(shoegaplabelloc)
            shoegaplabelloc(:,1) = -shoegaplabelloc(:,1);
        end
        
        if ~isempty(coillabellocs)
            coillabellocs(:,1) = -coillabellocs(:,1);
        end
        
        % rearrange the corner nodes to preserve clockwise ordering
        % starting from bottom left
        cornernodes = [cornernodes(2), cornernodes(1), cornernodes(4), cornernodes(3)];
        
    end
    
    % add the specified offset in the x direction
    nodes(:,1) = nodes(:,1) + xoffset;

    if ~isempty(shoegaplabelloc)
        shoegaplabelloc(:,1) = shoegaplabelloc(:,1) + xoffset;
    end

    if ~isempty(coillabellocs)
        coillabellocs(:,1) = coillabellocs(:,1) + xoffset;
    end
        
    elcount = elementcount_mfemm(FemmProblem);
    
    % store the nodes at the bottom of all the slots
    bottomnodes = cornernodes(1:2) + elcount.NNodes;  
    lastslottopnodes = cornernodes(3:4) + elcount.NNodes;    
    
    % draw the first slot linking it to the bottom of the domain
    thisslotsnodes = nodes;
        
    % move in the y direction to the first slot position
    thisslotsnodes(:,2) = thisslotsnodes(:,2) + slotpos(1);

    if ~isempty(coillabellocs)
        coillabellocs(:,2) = coillabellocs(:,2) + slotpos(1);
    end
        
    thisslotlinks = links + elcount.NNodes;
    
    [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, thisslotsnodes(:,1), thisslotsnodes(:,2));
        
    [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, thisslotlinks(:,1), thisslotlinks(:,2));
    
    if ~isempty(shoegaplabelloc)
        FemmProblem = addblocklabel_mfemm(FemmProblem, shoegaplabelloc(1,1), shoegaplabelloc(1,2) + slotpos(1), ...
                                            'BlockType', FemmProblem.Materials(Inputs.ShoeGapMaterial).Name, ...
                                            'MaxArea', Inputs.ShoeGapRegionMeshSize);
    end
    
    % draw the rest of the slots, making copies of the original slot's
    % nodes and links, adding in the labels, and linking them together
    for i = 2:numel(slotpos)
        
        thisslotsnodes = nodes;
        
        thisslotsnodes(:,2) = thisslotsnodes(:,2) + slotpos(i);
        
        thisslotlinks = links + numel(FemmProblem.Nodes);
        
        thisslotcornernodes = cornernodes + numel(FemmProblem.Nodes);
        
        [FemmProblem, nodeinds, nodeids] = addnodes_mfemm(FemmProblem, thisslotsnodes(:,1), thisslotsnodes(:,2));
        
        [FemmProblem, seginds] = addsegments_mfemm(FemmProblem, thisslotlinks(:,1), thisslotlinks(:,2));
        
        if ~isempty(shoegaplabelloc)
            
            FemmProblem = addblocklabel_mfemm(FemmProblem, shoegaplabelloc(1,1), shoegaplabelloc(1,2) + slotpos(i), ...
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
        coillabellocs = [coillabellocs; ...
                         coillabellocs(1:Inputs.NWindingLayers,1), ...
                         coillabellocs(1:Inputs.NWindingLayers,2) + slotpos(i) - slotpos(1)];
        
    end
    
    % we will return the outer corner node ids for later use
    outernodes = [bottomnodes, lastslottopnodes];
    
end


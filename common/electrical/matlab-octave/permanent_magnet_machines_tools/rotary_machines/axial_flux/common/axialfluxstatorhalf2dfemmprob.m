function [FemmProblem, outernodes, coillabellocs, slotinfo] = axialfluxstatorhalf2dfemmprob(slots, Poles, ypole, ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, xoffset, side, varargin)
% draw internal parts of half a slotted axial flux stator
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
    
    [FemmProblem, outernodes, coillabellocs, slotinfo] = ...
        uncurvedstatorhalf2dfemmprob ( numel(slotpos), ...
                                       ypole / slotsperpole, ...
                                       ycoil, ...
                                       yshoegap, ...
                                       xyoke, ...
                                       xcoil, ...
                                       xshoebase, ...
                                       xshoegap, ...
                                       xoffset, ...
                                       side, ...
                                      'NWindingLayers', Inputs.NWindingLayers, ...
                                      'FemmProblem', FemmProblem, ...
                                      'ShoeGapMaterial', Inputs.ShoeGapMaterial, ...
                                      'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                                      'Tol', Inputs.Tol, ... 
                                      'DrawCoilInsulation', Inputs.DrawCoilInsulation, ...
                                      'CoilInsulationThickness', Inputs.CoilInsulationThickness, ...
                                      'CoilBaseFraction', Inputs.CoilBaseFraction );
    
end


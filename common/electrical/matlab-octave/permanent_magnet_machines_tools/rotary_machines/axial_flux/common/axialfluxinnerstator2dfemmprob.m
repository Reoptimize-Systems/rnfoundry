function [FemmProblem, outernodes, coillabellocs] = axialfluxinnerstator2dfemmprob(yokecentresep, slots, Poles, ypole, ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, varargin)
% creates a femm problem description of an axial flux toothed armature
%
% Syntax
%
% [FemmProblem, outernodes, coillabellocs] = ...
%           axialfluxinnerstator2dfemmprob(yokecentresep, slots, Poles, ...
%                 ypole, ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, 'Parameter', Value)
%
% Input
%
%   yokecentresep - spacing between stator centres when drawing multiple
%     stators. Set to zero for a single stage machine.
%
%   slots - total number of slots in the machine
%
%   Poles - total number of Poles in the machine
%
%   ypole - pole pitch at the mean radius of the magnets
%
%   ycoil - coil pitch at the mean radius of the magnets
%
%   yshoegap - pitch of space between tooth shoes, i.e. the size of the
%     coil slot opening in the coil pitch direction.
%
%   xyoke - width of the stator yoke on which the slots are mounted
%
%   xcoil - width of the coil in the slot, i.e. the slot depth
%
%   xshoebase - width of the tooth shoe at the point where it joins the
%     tooth.
%
%   xshoegap - width of the tooth shoe at it's tip at the slot opening.
%
% A number of parameter-value paris can also be supplied (in the form of a
% string and the corresponding value) to set various optional arguments. If
% not supplied these are given default values. The possible parameters are
% listed below:
%
% 'FemmProblem' - an existing mfemm FemmProblem structure to which the new
%    elements will be added. If not supplied a new problem is created.
% 
% 'NStators' - integer number of stators to be drawn in total. Defaults to
%    one. If more than one are to be drawn they will drawn with the center
%    fo each yoke separated by the value of yokecentresep
%
% 'NWindingLayers' - integer number of axial coil layers in the design.
%
% 'SlotMaterial' - 
%
% 'SlotRegionMeshSize' - 
%
% 'ToothMaterial' - 
%
% 'ToothRegionMeshSize' - 
%
% 'ShoeGapMaterial' - 
%
% 'ShoeGapRegionMeshSize' - 
%
% 'NSlots' -  integer value determining the number of slots to be drawn. If
%   not supplied enough slots will be drawn to fill two Poles of the
%   machine design. This options can be used to draw large or smaller
%   simulations of the same design.
%
% 'Tol'
%
% Output
%
%   FemmProblem - An mfemm problem structure containing the new coil
%     elements.
%
%   outernodes - (NStators x 4) matrix containing the ids of the four outer
%     corner nodes of each drawn stator. The nodes are ordered in clockwise
%     direction starting from the bottom left corner node.
%
%   coillabellocs - (n * 2) matrix containing the x and y locations of the
%     coil labels. These are supplied for each stator part in sequence
%     moving from the left to right. 
%
%     [ 1st stator bottom slot outer coil layer x, 1st stator bottom left outer coil layer y ]
%     [ 1st stator bottom slot coil layer 2 x, 1st stator bottom left coil layer 2 y ]
%                                   .
%                                   . for NWindingLayers - 3 if more than 2
%                                   .
%     [ 1st stator bottom slot inner coil layer x, 1st stator bottom left inner coil layer y ]
%     [ 1st stator bottom left inner coil layer x, 1st stator bottom left inner coil layer y ]
%
% 

    Inputs.NStators = 1;
    Inputs.NWindingLayers = 1;
    Inputs.FemmProblem = newproblem_mfemm('planar');
    Inputs.SlotMaterial = 1;
    Inputs.SlotRegionMeshSize = -1;
    Inputs.ToothMaterial = 1;
    Inputs.ToothRegionMeshSize = -1;
    Inputs.ShoeGapMaterial = 1;
    Inputs.ShoeGapRegionMeshSize = -1;
    Inputs.NSlots = [];
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    Inputs = rmfield(Inputs, 'FemmProblem');
    
    stagepositions = (0:Inputs.NStators-1) .* yokecentresep;
    stagepositions = stagepositions - max(stagepositions)/2;
    
    coillabellocs = [];
    
    for i = 1:Inputs.NStators
        
        xoffset = stagepositions(i);
        
        % draw left side
        [FemmProblem, leftouternodes, leftcoillabellocs] = ...
            axialfluxstatorhalf2dfemmprob(slots, Poles, ypole, ycoil, ...
                      yshoegap, xyoke, xcoil, xshoebase, xshoegap, ...
                      xoffset, 'l', ...
                      'NWindingLayers', Inputs.NWindingLayers, ...
                      'FemmProblem', FemmProblem, ...
                      'ToothMaterial', Inputs.ToothMaterial, ...
                      'ToothRegionMeshSize', Inputs.ToothRegionMeshSize, ...
                      'ShoeGapMaterial', Inputs.ShoeGapMaterial, ...
                      'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                      'Tol', Inputs.Tol, ... 
                      'NSlots', Inputs.NSlots);
        
        % draw right side
        [FemmProblem, rightouternodes, rightcoillabellocs] = ...
            axialfluxstatorhalf2dfemmprob(slots, Poles, ypole, ycoil, ...
                      yshoegap, xyoke, xcoil, xshoebase, xshoegap, ...
                      xoffset, 'r', ...
                      'NWindingLayers', Inputs.NWindingLayers, ...
                      'FemmProblem', FemmProblem, ...
                      'ToothMaterial', Inputs.ToothMaterial, ...
                      'ToothRegionMeshSize', Inputs.ToothRegionMeshSize, ...
                      'ShoeGapMaterial', Inputs.ShoeGapMaterial, ...
                      'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                      'Tol', Inputs.Tol, ... 
                      'NSlots', Inputs.NSlots);
                  
                  
          outernodes(i,1:4) =  [leftouternodes(1), rightouternodes(2), rightouternodes(3), leftouternodes(4)]; 
              
          % store the coil label locations
          coillabellocs = [coillabellocs, [leftcoillabellocs; rightcoillabellocs]];
          
    end
    
end
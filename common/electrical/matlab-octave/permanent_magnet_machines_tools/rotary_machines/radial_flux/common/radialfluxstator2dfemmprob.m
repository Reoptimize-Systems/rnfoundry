function [FemmProblem, outernodes, coillabellocs] = radialfluxstator2dfemmprob(slots, poles, ryokecenter, thetapole, thetacoil, thetashoegap, ryoke, rcoil, rshoebase, rshoegap, drawnsides, varargin)
% creates a femm problem description of an axial flux toothed armature
%
% Syntax
%
% [FemmProblem, outernodes, coillabellocs] = ...
%           radialfluxinnerstator2dfemmprob(yokecentresep, slots, poles, ...
%                 thetapole, thetacoil, thetashoegap, ryoke, rcoil, rshoebase, rshoegap, 'Parameter', Value)
%
% Input
%
%   yokecentresep - spacing between stator centres when drawing multiple
%     stators. Set to zero for a single stage machine.
%
%   slots - total number of slots in the machine
%
%   poles - total number of poles in the machine
%
%   ryokecenter - radial position of the centre of the stator yoke (radial
%     distance from the centre of the machine).
% 
%   thetapole - pole pitch in radians 
%
%   thetacoil - coil pitch in radians
%
%   thetashoegap - pitch of space between tooth shoes, i.e. the size of the
%     coil slot opening in the coil pitch direction. Measure in radians.
%
%   ryoke - radial length of the stator yoke on which the slots are mounted
%
%   rcoil - radial length of the coil in the slot, i.e. the slot depth
%
%   rshoebase - radial length of the tooth shoe at the point where it joins
%     the tooth.
%
%   rshoegap - radial length of the tooth shoe at it's tip at the slot
%     opening.
%
%   drawnsides - 2 element vector determining which sides of a stator are
%     drawn. If the first element evaluates to true, an internally facing
%     stator is drawn, if the second element evaluates to true and
%     externally facing stator is drawn.
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
% 'NSlots' - integer value determining the number of slots to be drawn. If
%   not supplied enough slots will be drawn to fill two poles of the
%   machine design. This options can be used to draw large or smaller
%   simulations of the same design.
%
% 'Tol' - 
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

    if numel(drawnsides) ~= 2
        error('drawnsides must be a two element vector.')
    end
    
    Inputs.NWindingLayers = 1;
    Inputs.FemmProblem = newproblem_mfemm('planar');
    Inputs.SlotMaterial = 1;
    Inputs.SlotRegionMeshSize = -1;
    Inputs.ToothMaterial = 1;
    Inputs.ToothRegionMeshSize = -1;
    Inputs.ShoeGapMaterial = 1;
    Inputs.ShoeGapRegionMeshSize = -1;
    Inputs.ShoeGroup = 0;
    Inputs.YokeGroup = 0;
    Inputs.NSlots = [];
    Inputs.Tol = 1e-5;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = Inputs.FemmProblem;
    
    Inputs = rmfield(Inputs, 'FemmProblem');

    coillabellocs = [];
    
    if drawnsides(1)

        % draw inner internally facing side
        [FemmProblem, outercornernodes, outercoillabellocs] = ...
            radialfluxstatorhalf2dfemmprob(slots, poles, thetapole, thetacoil, ...
                      thetashoegap, ryoke, rcoil, rshoebase, rshoegap, ...
                      ryokecenter, 'i', ...
                      'NWindingLayers', Inputs.NWindingLayers, ...
                      'FemmProblem', FemmProblem, ...
                      'ToothMaterial', Inputs.ToothMaterial, ...
                      'ToothRegionMeshSize', Inputs.ToothRegionMeshSize, ...
                      'ShoeGapMaterial', Inputs.ShoeGapMaterial, ...
                      'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                      'ShoeGroup', Inputs.ShoeGroup, ...
                      'Tol', Inputs.Tol, ... 
                      'NSlots', Inputs.NSlots);
    end
    
    if drawnsides(2)
        
        % draw outer externally facing side
        [FemmProblem, innercornernodes, innercoillabellocs] = ...
            radialfluxstatorhalf2dfemmprob(slots, poles, thetapole, thetacoil, ...
                      thetashoegap, ryoke, rcoil, rshoebase, rshoegap, ...
                      ryokecenter, 'o', ...
                      'NWindingLayers', Inputs.NWindingLayers, ...
                      'FemmProblem', FemmProblem, ...
                      'ToothMaterial', Inputs.ToothMaterial, ...
                      'ToothRegionMeshSize', Inputs.ToothRegionMeshSize, ...
                      'ShoeGapMaterial', Inputs.ShoeGapMaterial, ...
                      'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                      'ShoeGroup', Inputs.ShoeGroup, ...
                      'Tol', Inputs.Tol, ... 
                      'NSlots', Inputs.NSlots);

    end
    
    if drawnsides(1) && drawnsides(2)
        outernodes =  [innercornernodes(1), outercornernodes(2), outercornernodes(3), innercornernodes(4)];
        coillabellocs = [innercoillabellocs; outercoillabellocs];
    elseif drawnsides(1)
        outernodes =  outercornernodes;
        coillabellocs = outercoillabellocs;
    elseif drawnsides(2)
        outernodes =  innercornernodes;
        coillabellocs = innercoillabellocs;
    else
        
    end
          
end
    

function FemmProblem = ax2rad(FemmProblem)

    for ind = 1:numel(FemmProblem.Nodes)
        [FemmProblem.Nodes(ind).Coords(1), FemmProblem.Nodes(ind).Coords(1)] = ...
            pol2cart(FemmProblem.Nodes(ind).Coords(2), FemmProblem.Nodes(ind).Coords(1));
    end
    
    for ind = 1:numel(FemmProblem.BlockLabels)
        FemmProblem.BlockLabels
    end

end
   

    
    
function [FemmProblem, outernodes, coillabellocs, inslabellocs] = radialfluxstatorhalf2dfemmprob(...
    slots, Poles, thetapole, thetacoil, thetashoegap, ryoke, rcoil, rshoebase, rshoegap, roffset, side, varargin)
% draw internal parts of a slotted stator for a radial flux machine
%
% Syntax
%
% [FemmProblem, outernodes, coillabellocs] = ...
%       radialfluxstatorhalf2dfemmprob( slots, Poles, thetapole, thetacoil, ...
%                                       thetashoegap, ryoke, rcoil, rshoebase, ...
%                                       rshoegap, coillayers, side )
%
%
% [FemmProblem, outernodes, coillabellocs] = ...
%       radialfluxstatorhalf2dfemmprob( ..., 'Parameter', Value )
%
% Description
%
% radialfluxstatorhalf2dfemmprob creates or adds to an mfemm FemmProblem
% structure, the internal parts of a radial flux slotted stator. The
% 'internal parts' refer to the segments, nodes and labels making up the
% teeth and coils, but not the outer part of the yoke of the machine, i.e.
% the stator drawing is not closed at its edges. It is intended for
% radialfluxstatorhalf2dfemmprob to be used by higher level drawing
% functions which are draing the full radial flux machine geometry. 
% 
% radialfluxstatorhalf2dfemmprob is a flexible function designed to draw
% any number of slots, but by default draws enough to create a drawing
% covering two full magnetic poles if possible.
%
% Input
%
%  slots - total number of slots in the machine design (note this is
%    different from the number of slots that will actually be drawn), used
%    to determine the slots / poles ratio.
%
%  Poles - total number of magnetic poles in the machine design (note this
%    is different from the number of poles that will actually be drawn),
%    used to determine the slots / poles ratio.
%
%  thetapole - 
%
%  thetacoil - 
%
%  thetashoegap - 
%
%  ryoke - 
%
%  rcoil - 
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
    Inputs.FemmProblem = newproblem_mfemm ('planar');
    Inputs.ShoeGapMaterial = 1;
    Inputs.ShoeGapRegionMeshSize = -1;
    Inputs.Tol = 1e-5;
    Inputs.CoilBaseFraction = 0.05;
    Inputs.ShoeCurveControlFrac = 0.5;
    Inputs.SplitX = false;
    Inputs.DrawCoilInsulation = false;
    Inputs.CoilInsulationThickness = 0;
    
    % check an integer number of machine Poles and slots (can't have half a
    % slot or a pole in a machine
    slots = round (slots);
    Poles = round (Poles);

    slotsperpole = slots / Poles;
    
    Inputs.NSlots = slotsperpole * 2;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    % copy over the FemmProblem
    FemmProblem = Inputs.FemmProblem;
    
    % remove the odl one to save memory
    Inputs = rmfield (Inputs, 'FemmProblem');
    
    [FemmProblem, outernodes, coillabellocs, slotinfo] = ...
        curvedstatorhalf2dfemmproblem ( Inputs.NSlots, ...
                                        thetapole / slotsperpole, ...
                                        thetacoil, ...
                                        thetashoegap,...
                                        ryoke, ...
                                        rcoil, ...
                                        rshoebase, ...
                                        rshoegap, ...
                                        roffset, ...
                                        side, ...
                                        'NWindingLayers', Inputs.NWindingLayers, ...
                                        'FemmProblem', FemmProblem, ...
                                        'ShoeGapMaterial', Inputs.ShoeGapMaterial, ...
                                        'ShoeGapRegionMeshSize', Inputs.ShoeGapRegionMeshSize, ...
                                        'Tol', Inputs.Tol, ... 
                                        'DrawCoilInsulation', Inputs.DrawCoilInsulation, ...
                                        'CoilInsulationThickness', Inputs.CoilInsulationThickness, ...
                                        'CoilBaseFraction', Inputs.CoilBaseFraction );
    
    inslabellocs = slotinfo.inslabelloc;
    
end


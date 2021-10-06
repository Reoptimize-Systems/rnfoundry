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
%   slots - total number of slots in the machine design (note this is
%    different from the number of slots that will actually be drawn), used
%    to determine the slots / poles ratio.
%
%   Poles - total number of magnetic poles in the machine design (note this
%    is different from the number of poles that will actually be drawn),
%    used to determine the slots / poles ratio.
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
%   roffset - radial position of the centre of the stator yoke (radial
%     distance from the centre of the machine).
%
%  In addition, a number of other optional arguments can be supplied as
%  parameter-value pairs. These options and their behaviour are as follows:
%
%   'FemmProblem' - an existing mfemm FemmProblem structure to which the new
%     elements will be added. If not supplied a new problem is created.
%
%   'NWindingLayers' - integer number of axial coil layers in the design.- 
%
%   'NSlots' - integer value determining the number of slots to be drawn. 
%     If not supplied enough slots will be drawn to fill two Poles of the
%     machine design. This options can be used to draw large or smaller
%     simulations of the same design.
%
%   'Tol' - tolerance at which to consider various dimensions to be zero,
%     by default this is 1e-5. This is used to prevent meshes occuring with
%     very large numbers of triangles.
%
%   'ShoeGapMaterial' - index of the material in the FemmProblem.Materials
%     structure containing the material to be used for the gap between teeth
%     shoes when the teeth have a shoe which ends in a blunt edge.
%
%   'ShoeGapRegionMeshSize' - Mesh size for the gap between teeth shoes when
%     the teeth have a shoe which ends in a blunt edge.
%
%   'CoilBaseFraction' - scalar value indicating the starting point of 
%     curvature at the base of the slot (closest to the yoke). By default
%     the coil slot is given a curved base. This value indicates where the
%     curvature begins, specified as a fraction of the distance between the
%     yoke and the start of the base of any shoe (if present), or just the
%     top of the slot, if no shoe is present. Default is 0.05 if not
%     supplied.
%
%   'ShoeCurveControlFrac' - factor controlling the 'curvature' of the 
%     tooth shoe, this is a value between 0 and 1. The exact effect of this
%     number is complex, and depends on the geometry of the slot. However,
%     in general a lower number results in a curve closer to a line draw
%     directly from the shoe base to the shoe gap, while higher numbers
%     aproximate a sharp right angle. Anything in between will produce a
%     smooth curve. Defaults to 0.5.
%
%     N.B. the slot geometry affects this curve in the following way. If
%     the position of the shoe gap node is below the intercept of the line
%     formed by the edge of the slot and a vertical line at the shoe gap
%     node, the resulting curve will bend outward from the inside of the
%     slot. If the intercept is below the shoe gap node, the curve will
%     bend into the slot.
%
%  'SplitSlot' - true/false flag. If there is only two winding layers, the 
%    slot can be split into two in the circumferential direction rather
%    than the radial by setting this flag to true. Defaults to false. If
%    true coil label locations are provided in an anti-clockwise direction.
%
%   'DrawCoilInsulation' = true/false flag indicating whether to draw a
%     layer of coil insulation in the slot
%
%   'CoilInsulationThickness' - scalar value giving the thickness of the
%     coil insulation to draw when DrawCoilInsulation is true. Default is 0
%     if not supplied, so no coil insulation will actually be drawn unless
%     you explicitly set a value greater than zero.
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
    Inputs.SplitSlot = false;
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
                                        'CoilBaseFraction', Inputs.CoilBaseFraction, ...
                                        'SplitSlot', Inputs.SplitSlot );
    
    inslabellocs = slotinfo.inslabelloc;
    
end


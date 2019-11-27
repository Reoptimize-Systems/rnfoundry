function total_joint = clampTwoNodes (node1, node2)
% clamp two nodes in initial relaitve position and orientation
%
% Syntax
%
%  total_joint = clampTwoNodes (node1, node2)
%
% Description
%
% clampTwoNodes is a convenience function which creates a total joint which
% clamps two nodes in their initial relative position and orientation.
%
% Input
%
%  node1 - mbdyn.pre.structuralNode (or derived class) object
%    representing to first node the joint connects
%
%  node2 - mbdyn.pre.structuralNode (or derived class) object
%    representing to second node the joint connects
%
% Output
%
%  total_joint - mbdyn.pre.totalJoint object with all position and 
%   orientation constraints active and set to keep the initial relative
%   position and orientation of the two input nodes.
%
%

    total_joint = mbdyn.pre.totalJoint ( node1, node2, ...
                                         'PositionStatus', [true, true, true], ...
                                         'OrientationStatus', [true, true, true], ...
                                         'ImposedRelativePosition', 'null', ...
                                         'RelativeOffset1', 'null', ...
                                         'RelativeOffset1Reference', 'node', ...
                                         'RelativeOffset2', 'null', ...
                                         'RelativeOffset2Reference', 'other node' );

end
function design = completedesign_FM(design, simoptions)
% completes common design aspects to create a minimal design of a flat profile
% linear machine
%
% Syntax
%
% design = completedesign_FM(design)
% design = completedesign_FM(design, simoptions)
%
% Description
%
% completedesign_FM adds various default design options common to all flat profile
% linear pm machines if they are not already present in the design structure. The
% following fields will be added to 'design' if they are not already present:
%
%  CoilLayers - default is 1
%
%  MagnetSkew - determines the amount of magnet skewwing as a ratio of a pole
%   width (i.e. it is expected to be between 0 and 1). Defaults to zero if not
%   supplied.
%
%  NStrands - number of strands making up the wire in the coils. Defaults to 1
%   if not supplied.
%
%  NStages - number of stages making up the machine. Defaults to 1 if not
%   supplied
% 
%  sides - provided for legacy code compatibility. Defaults to the same value as 
%   NStages if not supplied
%
% Not all machine types may make use of all the default options set here.
%
%
% See also: completedesign_linear.m, completedesign_AM.m
%

    if nargin < 2
        simoptions = struct ();
    end
    
    design = completedesign_linear(design, simoptions)
    
    design = setfieldifabsent(design, 'sides', design.NStages);
    
end
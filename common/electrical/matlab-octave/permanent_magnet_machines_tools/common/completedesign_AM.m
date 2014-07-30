function design = completedesign_AM(design, simoptions)
% completes common design aspects to create a minimal design ready for
% simulation
%
% Syntax
%
% design = completedesign_AM(design)
% design = completedesign_AM(design, simoptions)
%
% Description
%
% completedesign_AM adds various default design options common to all pm
% machines if they are not already present in the design structure. The
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
% Not all machine types may make use of all the default options set here.
%

    % set the magnet skew to zero by default
    design = setfieldifabsent (design, 'MagnetSkew', 0);

    % use single layered coils by default
    design = setfieldifabsent (design, 'CoilLayers', 1);
    
    % no coil insulation by default
    design = setfieldifabsent (design, 'CoilInsulationThickness', 0);
    
    % use single-stranded conductors by default
    design = setfieldifabsent (design, 'NStrands', 1);
    
    % check the number of stages in the design
    design = setfieldifabsent (design, 'NStages', 1);
    
end
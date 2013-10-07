function [vol, design] = structvol_ACPMSM(design, options)
% structvol_ACPMSM: calculates the total volume of the structure of the
% air-cored permanent magnet synchronous machine
%
% Includes the outer frame and guide rails in the calculation
%
% Syntax
%
% [vol, design] = structvol_ACPMSM(design, options)
%
    
    % set the FieldPoles parameter for the structural calculation
    design.FieldPoles = design.Poles(1);
    
    % calculate the depth of the machine for the calculation of web volume
    design.Depth = 2 * (design.dg + design.lm + design.dbi);
    
    % do the calculation using the common functions
    vol = structvol_FM(design, options);

end
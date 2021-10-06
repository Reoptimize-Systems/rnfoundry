function [vol, design] = structvol_PMSM(design, options)
% structvol_PMSM: calculates the total volume of the structure of the
% linear permanent magnet synchronous machine
%
% Includes the outer frame and guide rails in the calculation
%
% Syntax
%
% [vol, design] = structvol_PMSM(design, options)
%
    
    % set the FieldPoles parameter for the structural calculation
    design.FieldPoles = design.Poles(2);
    
    % calculate the depth of the machine for the calculation of web volume
    design.Depth = 2 * (design.hbf + design.hm + design.g + design.ht + design.hba);
    
    % do the calculation using the common functions
    vol = structvol_FM(design, options);
    
end
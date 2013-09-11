function vol = structvol_FM(design, options)
% structvol_FM: calculates the total volume of material used in the
% structure of a flat profile linear machine
%
% Includes the outer frame and guide rails in the calculation. It is always
% assumed it is a two-sided frame separated by web elements
%
% Syntax
%
% vol = structvol_FM(design, options)
%

    % First calculate the cross-sectional area of the beams
    area = CSArea(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod);
    
    % Next determine the volume of a single support beam
    volume = area .* design.ls .* options.alphab;
    
    n_beams = numbeams(span1(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod), ...
                       design.PoleWidth, ...
                       design.FieldPoles * design.PoleWidth, ...
                       design.BeamSpreadFactor);
    
    % we always assume there are two sides to the frame
    vol = 2 .* volume .* n_beams;
    
    % Now add the supporting beams to which the pole supports are fixed,
    % these are 
    height = (design.FieldPoles .* design.PoleWidth) - design.PoleWidth;
    
    vol = vol + (4 * CSArea(design.OuterStructureBeamVars(2,:), ...
                            design.OuterStructureBeamIMethod) * height);
                        
    % add the outer webs holding the two sides of the frame apart
    vol = vol + (design.OuterWebs * CSArea(design.OuterStructureBeamVars(2,:), ...
                                           design.OuterStructureBeamIMethod) ...
                                   * design.Depth);
                                     
                        
	% finally, add the guide rails to which the frame is fixed
	height = max(design.poles) .* design.PoleWidth;
    
	vol = vol + (4 * CSArea(design.GuideRailIVars, ...
                            design.GuideRailIMethod) * height);

end
function BeamInfo = completebeaminfo_FM(design, options, BeamInfo)
    
    %Young's modulus of the beam material in Pa
    BeamInfo.SupportBeams.E = options.E(1);

    % Poisson ratio
    BeamInfo.SupportBeams.nu = 0.31;

    % Get the properties of the beam section
    BeamInfo.SupportBeams.A = CSArea(design.OuterStructureBeamVars(2,:), design.OuterStructureBeamIMethod);
    BeamInfo.SupportBeams.I2 = MomentOfInertiaY1(design.OuterStructureBeamVars(2,:), design.OuterStructureBeamIMethod);
    BeamInfo.SupportBeams.I3 = MomentOfInertiaY2(design.OuterStructureBeamVars(2,:), design.OuterStructureBeamIMethod);
    BeamInfo.SupportBeams.I1 = BeamInfo.SupportBeams.I2 + BeamInfo.SupportBeams.I3;
    % Temporarily set 'J' the torsion constant to be the same as the polar
    % moment of inertia. This will result in an overestimate of the
    % stiffness in torsion
    BeamInfo.SupportBeams.J = BeamInfo.SupportBeams.I1;
    BeamInfo.SupportBeams.rho = options.StructMaterialDensity;
    BeamInfo.SupportBeams.AngleFromHorizontal = design.AngleFromHorizontal;
    BeamInfo.SupportBeams.NoPerSide = 2; 

    % Initialise the guide rails info with values from the support beams,
    % there will be some common properties
    BeamInfo.GuideRails = BeamInfo.SupportBeams;
    BeamInfo.GuideRails.A = CSArea(design.GuideRailIVars, design.GuideRailIMethod);
    BeamInfo.GuideRails.I2 = MomentOfInertiaY1(design.GuideRailIVars, design.GuideRailIMethod);
    BeamInfo.GuideRails.I3 = MomentOfInertiaY2(design.GuideRailIVars, design.GuideRailIMethod);
    BeamInfo.GuideRails.I1 = BeamInfo.GuideRails.I2 + BeamInfo.SupportBeams.I3;
    BeamInfo.GuideRails.J = BeamInfo.GuideRails.I1;
    BeamInfo.GuideRails.length = max((design.NInnerPoles .* design.PoleWidth) - design.PoleWidth, design.PoleWidth);
    BeamInfo.GuideRails.Sections = 5;
    
    % Set up the bearing connections between the guide rails and the frame
    BeamInfo.GuideBearings = BeamInfo.SupportBeams;
    BeamInfo.GuideBearings.NoPerGuide = 2;
    % Make the guide rail->field connections super stiff
    BeamInfo.GuideBearings.E = BeamInfo.GuideBearings.E * 100;

    % Set up the support beams to which the outer pole supports are
    % attached
    BeamInfo.OuterPoleSupports = BeamInfo.SupportBeams;
    BeamInfo.OuterPoleSupports.A = CSArea(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod);
    BeamInfo.OuterPoleSupports.I2 = MomentOfInertiaY1(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod);
    BeamInfo.OuterPoleSupports.I3 = MomentOfInertiaY2(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod);
    BeamInfo.OuterPoleSupports.I1 = BeamInfo.OuterPoleSupports.I2 + BeamInfo.OuterPoleSupports.I3;
    BeamInfo.OuterPoleSupports.J = BeamInfo.OuterPoleSupports.I1;
    BeamInfo.OuterPoleSupports.rho = options.StructMaterialDensity;
    BeamInfo.OuterPoleSupports.Sections = max(4, round2even(design.ls * options.sections));
    BeamInfo.OuterPoleSupports.NoPerSide = max(2, numbeams(span1(design.OuterStructureBeamVars(1,:), design.OuterStructureBeamIMethod), ...
        design.PoleWidth, design.NOuterPoles .* design.PoleWidth, design.BeamSpreadFactor));
    
    % Set up the outer sturcture webs, i.e. the structural members which
    % keep the two parts of the outer structure apart
    BeamInfo.OuterWebs = BeamInfo.SupportBeams;
    BeamInfo.OuterWebs.NoPerSide = design.OuterWebs;
    BeamInfo.OuterWebs.Sections = max(3, round2even(BeamInfo.depth * options.sections));
    
    BeamInfo.PreviousMeshes = {};
    BeamInfo.PreviousMeshesNBeams = [];
    
end
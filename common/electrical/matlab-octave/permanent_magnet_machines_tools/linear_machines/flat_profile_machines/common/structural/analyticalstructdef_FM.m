function [Def, BeamInfo] = analyticalstructdef_FM(design, options, WperRail, BeamInfo, forceRatio)
            
    % get the starting air gap closing forces 
    
    % The following assumes a polynomial has been fitted to
    % the force per unit area of translator surface versus
    % the air gap
    closingForce = -polyvaln(design.p_gforce, design.g);

    % Get the force per unit length on a support beam
    closingForce = closingForce * design.PoleWidth * forceRatio;
    
    % Test deflection of beam between support points
    % Outer pole supports rail will be modelled as beam fixed at both
    % ends with a superposition of concentrated loadings
    
    % find longest span between supports of any kind on outer frame
    [x, y, z, n, ...
        znodes, ...
        zguidecon, ...
        zsupp, ...
        zframe, ...
        zextreme, ...
        guideyshift, ...
        tolerance] = dimsfebeamdef_FM(BeamInfo);

    temp = uniquetol([zsupp; zguidecon], tolerance);
    
    temp = sortrows(temp);
    
    BeamInfo.SupportBeams.MaxSpan = max(temp(2:end) - temp(1:end-1));
    
    supportsInSpan = ceil((BeamInfo.SupportBeams.MaxSpan / BeamInfo.height) * BeamInfo.OuterPoleSupports.NoPerSide);
    
    supportPositions = (0:(BeamInfo.SupportBeams.MaxSpan / (supportsInSpan + 1)):BeamInfo.SupportBeams.MaxSpan)';
 
    beamMethod = '3.1d';

    Yvars = zeros(numel(supportPositions)-2, 3);

    Yvars(:,1) = closingForce * design.ls * options.alphab / 2;

    Yvars(:,2) = BeamInfo.SupportBeams.MaxSpan;

    Yvars(:,3) = supportPositions(2:end-1);

    defpoints = linspace(0, Yvars(1,2), 20);

    Def = max(BeamDeflectionSuperY1(Yvars, design.OuterStructureBeamVars(2,:), BeamInfo.GuideRails.E, defpoints, design.OuterStructureBeamIMethod, beamMethod));

    % Guide rail will be modelled as beam simply supported at both ends
    beamMethod = '3.1d';
    
    %   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - W
%          Yvars(:,2) - l, length of the beam
%          Yvars(:,3) - a, distance from M_A at which 'W' is applied 

    Yvars = zeros(numel(zguidecon), 3);

    Yvars(:,1) = WperRail / numel(zguidecon);

    Yvars(:,2) = 2 * zextreme;

    Yvars(:,3) = zextreme + zguidecon;

    defpoints = linspace(0, Yvars(1,2), 20);

    Def = Def + max(BeamDeflectionSuperY1(Yvars, design.GuideRailIVars, BeamInfo.GuideRails.E, defpoints, design.GuideRailIMethod, beamMethod));

    % Outer pole support beams will be modelled as beam fixed at both
    % ends
    beamMethod = '3.2d';

    Yvars = zeros(1, 4);

    Yvars(:,1) = closingForce;

    Yvars(:,2) = closingForce;

    Yvars(:,3) = design.ls * options.alphab;

    Yvars(:,4) = 0;

    defpoints = linspace(0, Yvars(1,3), 20);

    Def = Def + max(BeamDeflectionSuperY1(Yvars, design.OuterStructureBeamVars(1,:), BeamInfo.GuideRails.E, defpoints, design.OuterStructureBeamIMethod, beamMethod));

    
end
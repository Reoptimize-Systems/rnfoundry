function [design, PorF, z, peakI, peakJz, Def, fFVars, aFVars, suppWperm, polesPerSection, a1, a2, fstart, fend] = evaluatestructure_TM(design, options)

    PorF = [0 0 0 0];
    
    % 0. Electrical info such as peak current etc will be passed into the
    % function

    % 1. Determine forces in coil based on max allowable current, e.g. more
    % than current density of 5e6 A/m^2. The maximum force occurs 

    peakI = 5e6 * (design.Dc/2)^2 * pi;
    peakJz = peakI .* design.CoilTurns ./ design.CoilArea;
    
%     [Force, avFro, avFri] = coilgapclosingforce_TM(design, peakJz, (design.Wp/8) ./ design.Wp, 0, 1);
    [avFro, avFri] = avRoandRiforces(design, peakJz);
    
    % 2. Determine max stresses in coil
    %   cVars - (n x 3) matrix:-
    %           Col 1. a, the outer radius
    %           Col 2. b, the inner radius
    %           Col 3. v, Poisson's ratio for the material  
    cVars = [design.Ro, design.Rm, options.vCoil];
    MaxSigma2 = maxcylstressradialforce(avFri, avFro, cVars);

    % 3. Does it exceed the yield point?
    % Return the maximum coil stresses, a higher function will determine
    % the correct course of action based on this value
    PorF(1) = abs(MaxSigma2);
    
    % 4. Determine lateral forces in offset coil worst case scenario, these
    % will be zero though!! Perhaps use wave loading force or own weight at 45
    % degree incline or something    
    
    % Get the initial deflections due to manufacturing tolerance in the
    % sheath and outer diameter of the shaft  
    
%     sheathThickness = design.Ra - design.Ro;
    
    % 10% or 0.2 mm tolerance on the thickness of the tube wall based on
    % british standard BS EN 10297-2 for seamless cold-rolled stainless
    % steel tubes. I saved this in a word file somewhere
    
%     RatTol = sheathThickness * 0.1;
%     
%     if RatTol < 0.2/1000 
%         RatTol = 0.2/1000;
%     end
    RatTol = 0;
    
%     % 0.75 % or 0.3 mm tolerance on the outer diameter 
%     RaoDTol = design.Ra * 0.0075;
%     if RaoDTol < 0.3/1000
%        RaoDTol =  0.3/1000;
%     end

%     % 0.75 % or 0.3 mm tolerance on the inner diameter (whatever tube the 
%     % coil's wound on) 
    RaoDTol = design.Ri * 0.0075;
%     if RaoDTol < 0.3/1000
%        RaoDTol =  0.3/1000;
%     end
%     RaoDTol = 0;
% 
%     soDTol = design.Rso * 0.0075;
%     if soDTol < 0.3/1000
%         soDTol = 0.3/1000;
%     end
    soDTol = 0;
%     
%     RmoDTol = design.Rm * 0.0075;
%     if RmoDTol < 0.3/1000
%        RmoDTol =  0.3/1000;
%     end
%     RmoDTol = 0;

    % 0.2 mm standard tolerance on outer diameter of steel and magnet discs
    RmoDTol = 0.2 / 1000;
    
    % determine locations of points between sections, choosing an even
    % number of sections will ensure the centre of the beam is chosen
    
    % the sections will all be in the active part where most loading will
    % occur, if there are supports we will have 3 sections at either end of
    % the translator where only the supports are present (the armature
    % supports being identical to the rest of the armature)
    a2 = design.totalLength(1) - design.supportLengths(1,2);
    % Next determine the position of the start of the active part
    a1 = design.supportLengths(1,1);
    
    % Make some deflection sample points along the first support 
    if a1 > 0
        supSecs = [3, 0];
        z = 0:a1/supSecs(1):a1;
        fstart = supSecs(1) + 1;
    else
        supSecs = [0, 0];
        z = 0;
        fstart = 1;
    end
    
    fend = supSecs(1) + options.StructSections;
    
    fieldLen = design.totalLength(1) - sum(design.supportLengths(1,:));
    
    % Add in the points along the active part
    z = [z, (a1+fieldLen/options.StructSections):(fieldLen/options.StructSections):a2];

    % Add in the second set of support points
    if a2 >= design.totalLength(1)
%         z = [z, (a2+(design.totalLength(1)-a2)/supSecs):(design.totalLength(1)-a2)/supSecs:design.totalLength(1)];
        supSecs(2) = 0;
    else
        supSecs(2) = 3;
        z = [z, (a2+(design.totalLength(1)-a2)/supSecs(2)):(design.totalLength(1)-a2)/supSecs(2):design.totalLength(1)];
%         z = [z, a2:(design.totalLength(1)-a2)/supSecs(2):design.totalLength(1)];
    end
    
    % Initialise the deflections at these points to the tolerances defined
    % previously
    Def = [repmat((RaoDTol+RatTol), 1, size(z,2)); repmat(soDTol + RmoDTol, 1, size(z,2))];

    % Initialise forces to own weight at given incline, in this case we
    % will consider the field forces to be negative. This is because the
    % deflections are always calculated as positive since we are summing
    % them to get the total deflection. In this case though, both field and
    % armature will be bending in the same direction, so the deflections
    % will cancel out rather than reinforce
    if ~isfield(design, 'AngleFromHorizontal')
        design.AngleFromHorizontal = pi/4;
    end

    if all(isfield(design, {'Rs2VHmag', 'Rs1VHmag', 'Ws2VhalfWs', 'Ws1VhalfWs'}))
        [design.fpW] = fieldpoleweight_TM(design.WmVWp, design.WpVRm, ...
                                      design.RsiVRso, design.RsoVRm, ...
                                      design.Rm, options.FieldIronDensity, ...
                                      options.MagnetDensity, options.StructMaterialDensity, ...
                                      design.Rs2VHmag, design.Rs1VHmag, ...
                                      design.Ws2VhalfWs, design.Ws1VhalfWs) * cos(design.AngleFromHorizontal);
    else
        [design.fpW] = fieldpoleweight_TM(design.WmVWp, design.WpVRm, ...
                                      design.RsiVRso, design.RsoVRm, ...
                                      design.Rm, options.FieldIronDensity, ...
                                      options.MagnetDensity, ...
                                      options.StructMaterialDensity) * cos(design.AngleFromHorizontal);
    end
    
    [design.apW] = armaturepoleweight_TM(design.WpVRm, design.RoVRm, design.Rm, ...
                                  design.g, design.WcVWp, design.fillfactor, ...
                                  options.CopperDensity, design.RaVRo, ...
                                  options.ArmatureIronDensity) * cos(design.AngleFromHorizontal);
    
    fFVars = repmat((-design.fpW./design.Wp)*(fieldLen/options.StructSections),1,options.StructSections);% zeros(size(x,2));
    aFVars = repmat((design.apW/design.Wp)*(fieldLen/options.StructSections),1,options.StructSections);
    
    % Calculate the weight per m of the support parts, in the case of the
    % field this is the weight per m of the shaft, but in the case of the
    % armature this is the weight per m of the combined coil and sheath.
    % Note that the armature pole weight has already been adjusted to
    % account fo the design.AngleFromHorizontal, so we do not need to do this again, but we must
    % do this for the shaft support, again the field force is negative
    suppWperm = [-(pi * (design.Rso^2 - design.Rsi^2) * options.StructMaterialDensity * cos(design.AngleFromHorizontal) * 9.81), (design.apW / design.Wp)] ;
    
    polesPerSection = fieldLen / (design.Wp * options.StructSections);
    
end

function [avFro, avFri] = avRoandRiforces(design, Jz)

    y = (design.Wp/8) ./ design.Wp;
    ymin = y - (design.WcVWp / 2);
    ymax = y + (design.WcVWp / 2);

    avFro = quad(@(z) lineBy(design.Ro, z, design.p_By), ymin, ymax) * design.Wp .* Jz  .* (pi * 2 * design.Ro) ./ (pi * 2 * design.Ro * design.Wc );

    avFri = quad(@(z) lineBy(design.Ri, z, design.p_By), ymin, ymax) * design.Wp .* Jz .* (pi * 2 * design.Ri) ./ (pi * 2 * design.Ri * design.Wc );

end

function By = lineBy(xtemp, ytemp, p_By)

    xtemp = repmat(xtemp, size(ytemp));

    By = polyvaln(p_By, [xtemp', ytemp'])';

end

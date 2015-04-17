function FemmProblem = singleslot_RADIAL_SLOTTED (design)
% draws nodes arcs and segments for a single slot in a slotted radial flux
% machine
%
% Syntax
%
% FemmProblem = singleslot_RADIAL_SLOTTED (design)
%

    switch design.ArmatureType
        
        case 'external'
            drawnstatorsides = [1, 0];
            Rs = design.Rmo + design.g + design.tc(1) + design.tsb + design.ty/2;
        case 'internal'
            drawnstatorsides = [0, 1];
            Rs = design.Rmi - design.g - design.tc(1) - design.tsb - design.ty/2;
            
    end
    
    if numel (design.tc) > 1
        coilbasefrac = design.tc(2) / design.tc(1);
    else
        coilbasefrac = 0.05;
    end
    
    if drawnstatorsides(1)

        % draw inner internally facing side
        FemmProblem = ...
            radialfluxstatorhalf2dfemmprob(design.Qs, design.Poles, design.thetap, design.thetac, ...
                      design.thetasg, design.ty, design.tc(1), design.tsb, design.tsg, ...
                      Rs, 'i', ...
                      'NWindingLayers', design.CoilLayers, ...
                      'NSlots', 1, ...
                      'CoilBaseFraction', coilbasefrac);
    end
    
    if drawnstatorsides(2)
        
        % draw outer externally facing side
        FemmProblem = ...
            radialfluxstatorhalf2dfemmprob(design.Qs, design.Poles, design.thetap, design.thetac, ...
                      design.thetasg, design.ty, design.tc(1), design.tsb, design.tsg, ...
                      Rs, 'o', ...
                      'NWindingLayers', design.CoilLayers, ...
                      'NSlots', 1, ...
                      'CoilBaseFraction', coilbasefrac);

    end
    
end
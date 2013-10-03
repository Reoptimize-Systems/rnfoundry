function design = designstructure_ACTIAM(design, options, Sr, Sz)
% designstructure_ACTIAM: chooses the optimum inner shaft radius for a
% machine to minimize cost and generates design.StructPenalty factors where the design
% is not feasible 
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RaVRo - scalar value of Ra/Ro Ratio, the sheath outer radius to inner
%           radius ratio
%
%   RsoVRm - scalar value of the outer shaft radius to the translator
%            radius, Rso / Rm
%
%   WcVWp - scalar ratio of the coil width to pole width in the machine 
%           (Wc / Wp)
%
%   Rm - Translator radius in m
%
%   g - the machine aripag size
%
%   Ntot - total number of turns in coil
%  
%   dc - the conductor diameter in the winding
%
%   supportLengths - (2 X 2) matrix of lengths of the supports at either
%                    end of the active part of the field/translator. I.e.
%                    the first value is the distance from the first bearing
%                    to the start of the active part of the translator,
%                    while the second is the distance from the end of the
%                    active part to the second bearing.
%
%   totalLength - (1 X 2) Vector. The total length of the field and
%                 armature support respectively.
%
%   sections - 
%
%   Sr - the number of sections the coil is split into for evaluation
%        radially
%
%   Sz - the number of sections the coil is split into for evaluation
%        axially
%
% Output:
%
%   design.StructPenalty - vector containig two values, the first is any design.StructPenalty
%             assciated with coil failure in hoop. The second is any
%             design.StructPenalty associated with beam failure in both the field and
%             armature
%
%   RsiVRm - the optimal value of RsiVRm to achieve the given
%            specifications
%
%   design.actForce - the total force acting between the points in z at the point
%              at which the algorithm terminated
%

% Copyright Richard Crozer, The University of Edinburgh

    if ~isfield(design, 'defThresh')

        design.defThresh = 0.1 * design.g;

    end
    
    design.StructPenalty = [0 0];

    % First test the armature with a super stiff solid shaft to see if
    % armature fails as there's not much we can do in this case
    design.RsiVRso = 0;

    E = [207e9 100e6 207e9];

    Etemp = E;
    Etemp(1) = E(3) * 10000;

    [PorF, z, design.actForce, design.maxDefCriteria, netResultantForce, P] = ...
        evaluatestructure_ACTIAM(design, options, Etemp);
    design.StructDeflections = PorF;
        
    % if the coil hoop stress exceeds the maximum acceptable, generate a
    % penaty factor for the design
    if PorF(1) > options.coilYieldStrength
        design.StructPenalty(1) = (PorF(1)-options.coilYieldStrength)/options.coilYieldStrength;
    end

    % If the armature fails then there is no configuration which can
    % succeed so we add a design.StructPenalty factor. We will also return
    % a solid shaft which will add extra cost, but we will test this first
    % to see if this would also fail and add extra penalties if this is the
    % case
    if PorF(4) > design.defThresh
        design.StructPenalty(2) = (PorF(4)-design.defThresh) ./ (design.g - design.defThresh);
    end
    
    % Next test the shaft with solid cross section while making the
    % armature super stiff relative to it to see if there's any possibility
    % that any shaft can work
    Etemp = E;
    Etemp(3) = E(1) * 10000;
    
    [PorF, z, design.actForce, design.maxDefCriteria] = evaluatestructure_ACTIAM(design, options, Etemp);
    design.StructDeflections = [design.StructDeflections; PorF];
    
    % If the shaft also fails on it's own, add a further
    % design.StructPenalty to account for this
    if PorF(4) > design.defThresh
        design.StructPenalty(2) = design.StructPenalty(2) + (PorF(4)-design.defThresh) ./ (design.g - design.defThresh);
    end
    
    % if either part has totally failed return the current values
    if design.StructPenalty(2) > 0
        return;
    end
    
    % If neither part has failed individually we test the normal strength
    % armature and solid shaft together to see if the gap is maintained
    % with the strongest shaft available
    [PorF, z, design.actForce, design.maxDefCriteria, netResultantForce, P] = ...
        evaluatestructure_ACTIAM(design, options, E, netResultantForce, P);
    design.StructDeflections = [design.StructDeflections; PorF];
    
    % if the threshold is breached with the most solid shaft available when
    % both components are allowed to have the correct stiffness we add a
    % design.StructPenalty factor and return the solid shaft                                             
    if PorF(4) > design.defThresh
        % In this case the design.StructPenalty factor is the factor by
        % which the total deflection exceeds the allowed deflection
        design.StructPenalty(2) = (PorF(4)-design.defThresh) ./ (design.g - design.defThresh);
    end
    
    % if the parts have failed together return the current values
    if design.StructPenalty(2) > 0
        return;
    end
    
    % If we are stil in business, attempt to find the minimum shaft size
    % that will succeed. We will use fminbound to find this, as we know the
    % function will be fairly well behaved due to the previous testing.
    x1 = 0; 
    x2 = 0.95;
    % This will find the RsiVRm value that will give the closest deflection
    % to the deflection threshold
    design.RsiVRso = fminbnd(@defThreshDiff,x1,x2);
    design = ratios2dimensions_ACTIAM(design);
    design.StructDeflections = [design.StructDeflections; PorF];
    
    function diff = defThreshDiff(tempRsiVRso)
        
        design.RsiVRso = tempRsiVRso;
        design = ratios2dimensions_ACTIAM(design);
        
        [PorF, z, design.actForce, design.maxDefCriteria, netResultantForce, P] = ...
            evaluatestructure_ACTIAM(design, options, E, netResultantForce, P);
        
        diff = abs(design.defThresh-abs(PorF(4)));
        
    end
    
                                           
end
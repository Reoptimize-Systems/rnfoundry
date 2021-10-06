function varargout = evaluatestructure_ACPMSM(IVars, IMethod, design, E, sections, initDef, alphab, gfactor)
% evaluatestructure_ACPMSM: evaluates whether a given beam can withstand
% the forces in an air-cored linear permanent magnet synchronous machine
% design
%
% 
    design = ratios2dimensions_ACPMSM(design);
                                                  
    % We will calculatte the deflection of the beams, stopping when we reach
    % an beam that meets the specifications, as this will be the smallest
    % volume with which this can be achieved
              
    if strcmp(IMethod, '1.3')
        forceRatio = forceratio_linear(design, IVars(2));
    elseif strcmp(IMethod, '1.6')
        forceRatio = forceratio_linear(design, IVars(2));
    else
        error('Only rectangular and I-Beam sections are currently supported')
    end
    
    % Determine what the air-gap is with the structure loaded by the maxwell
    % stresses
    new_g = airgapclosure_ACPMSM(design, E, IVars, sections, IMethod, forceRatio, initDef, alphab);

    % Check if beam is successful or a failure
    if min(min(new_g)) >= gfactor * design.g
        varargout{1} = 1;
        varargout{2} = IVars;
    else
        varargout{1} = 0;
        varargout{2} = IVars;
    end

end
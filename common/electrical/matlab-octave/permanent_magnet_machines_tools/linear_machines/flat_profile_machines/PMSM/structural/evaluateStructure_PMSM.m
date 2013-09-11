function varargout = evaluateStructure_PMSM(IVars, IMethod, design, E, sections, mode, initDef, alphab, gfactor)
% evaluateStructure_PMSM: evaluates whether a given beam can withstand
% the forces in a linear permanent magnet synchronous machine design
%
% 
    design = ratios2dimensions_PMSM(design);
                                                  
    % We will calculatte the deflection of the I-Beams, stopping when we reach
    % an I-Beam that meets the specifications, as this will be the smallest
    % volume with which this can be achieved

%   IVarsCell - (1 x 2) cell array of values, each of which contains a matrix
%           of values for calculating the second moment of area of the
%           stator and translator sections repectively, the the values in
%           the second cell are used for both stator sections if a double
%           sided machine is being investigated.
%
%           First cell of IVars must be the beams supporting the inner and
%           outer 
%
%           The second cell of IVars
%
%           IVars{2}(1,1): dbi, half thickness of central section (back 
%                          iron thickness)
%           IVars{2}(1,2): Taup, the pole width
%           IVars{2}(1,3): hs, the tooth height
%           IVars{2}(1,4): bt, the tooth width
    
    % complete the IVars cell array by adding the appropriate information
    % from the design structure
    IVarsCell = {IVars, [design.hba,  design.Wp,   design.ht,    design.Wt]};
              
    if strcmp(IMethod, '1.3')
        forceRatio = forceratio_linear(design, IVars(2));
    elseif strcmp(IMethod, '1.6')
        forceRatio = forceratio_linear(design, IVars(2));
    else
        error('Only rectangular and I-Beam sections are currently supported')
    end
    
    % Determine what the air-gap is with the structure loaded by the maxwell
    % stresses
    new_g = airgapclosure_PMSM(design, E, IVarsCell, IMethod, sections, mode, forceRatio, initDef, alphab);

    % Check if beam is successful or a failure
    if min(min(new_g)) >= gfactor * design.g
        varargout{1} = 1;
        varargout{2} = IVars;
    else
        varargout{1} = 0;
        varargout{2} = IVars;
    end

end
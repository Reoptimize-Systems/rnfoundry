function [design, simoptions] = evaluatestructure_RADIAL_SLOTTED(design, simoptions)
% evaluates the structure of a slotted radial flux machine
%
% Syntax
%
% [design, simoptions] = evaluatestructure_RADIAL_SLOTTED(design, simoptions)
%
% 

    if strcmp(design.ArmatureType, 'internal')
        Iomega = 0;
        Oomega = design.OmegaPeak;
    elseif strcmp(design.ArmatureType, 'external')
        Iomega = 0;
        Oomega = design.OmegaPeak;
    else
        error('Only ''internal'' and ''external'' stator types supported.')
    end
    
    % evaluate the inner structure of the machine
    [mingapI, maxstressI, design] = ...
        evaluateinnerstructure_RADIAL_SLOTTED(design, simoptions, Iomega);
    
    % evaluate the outer strucutre of the machine
    [mingapO, maxstressO, design] = ...
        evaluateouterstructure_RADIAL_SLOTTED(design, simoptions, Oomega);
    
    % find the minimum air gap due to both inner and outer structure
    % deflections
    design.Min_g = design.g - ((design.g - mingapI) + (design.g - mingapO));
    
    % find the maximum stress in either the outer or inner structure
    design.Max_Stress = max(maxstressO, maxstressI);
    
end
function [design, simoptions] = structuralpenalty_FM(design, simoptions)

    % exceeding max allowed deflection
    design.maxAllowedDeflectionpenalty = 0;
    if design.StructPenalty(1) > 0
        % default to a linear penalty
        simoptions = setfieldifabsent(simoptions, 'maxAllowedDeflectionPenFactor', [1, 0]);

        if isscalar(simoptions.maxAllowedDeflectionPenFactor)
            simoptions.maxAllowedDeflectionPenFactor = [simoptions.maxAllowedDeflectionPenFactor, 0];
        end

        design.maxAllowedDeflectionpenalty = ...
            simoptions.maxAllowedDeflectionPenFactor(1) * design.StructPenalty(1)  ...
            +  (simoptions.maxAllowedDeflectionPenFactor(2) * design.StructPenalty(1))^2;

    end
        
end
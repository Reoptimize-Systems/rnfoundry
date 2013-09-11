function [score, design] = masspenalties_ROTARY(design, simoptions, score)

    if nargin < 3
        score = 0;
    end
    
    % exceeding max allowed rotor mass
    design.maxAllowedRotorMassPenalty = 0;

    if isfield(simoptions, 'maxAllowedRotorMass')
        if ~isempty(simoptions.maxAllowedRotorMass)
            if design.RotorMass > simoptions.maxAllowedRotorMass
                
                design.maxAllowedRotorMassPenalty = ...
                    10 * design.BaseScore * (design.RotorMass/simoptions.maxAllowedRotorMass) + ...
                     1 * design.BaseScore * (10 * design.RotorMass/simoptions.maxAllowedRotorMass)^2;

                score = score + design.maxAllowedRotorMassPenalty;
                
            end
        end
    end

    design.maxAllowedRotorModuleMassPenalty = 0;

    % exceeding max allowed rotor module mass
    if isfield(simoptions, 'maxAllowedRotorModuleMass')
        if ~isempty(simoptions.maxAllowedRotorModuleMass)
            if design.RotorModuleMass > simoptions.maxAllowedRotorModuleMass
                
                design.maxAllowedRotorModuleMassPenalty = ...
                    10 * design.BaseScore * (design.RotorModuleMass/simoptions.maxAllowedRotorModuleMass) + ...
                     1 * design.BaseScore * (10 * design.RotorModuleMass/simoptions.maxAllowedRotorModuleMass)^2;

                score = score + design.maxAllowedRotorModuleMassPenalty;
                
            end
        end
    end

    % exceeding max stator mass
    design.maxAllowedStatorMassPenalty = 0;

    if isfield(simoptions, 'maxAllowedStatorMass')
        if ~isempty(simoptions.maxAllowedStatorMass)
            if design.StatorMass > simoptions.maxAllowedStatorMass
                
                design.maxAllowedStatorMassPenalty = ...
                    10 * design.BaseScore * (design.StatorMass/simoptions.maxAllowedStatorMass) + ...
                     1 * design.BaseScore * (10 * design.StatorMass/simoptions.maxAllowedStatorMass)^2;

                score = score + design.maxAllowedStatorMassPenalty;
                
            end
        end
    end

    design.maxAllowedStatorModuleMassPenalty = 0;

    % exceeding max stator module mass
    if isfield(simoptions, 'maxAllowedStatorModuleMass')
        if ~isempty(simoptions.maxAllowedStatorModuleMass)
            if design.RotorModuleMass > simoptions.maxAllowedStatorModuleMass
                
                design.maxAllowedStatorModuleMassPenalty = ...
                    10 * design.BaseScore * (design.StatorModuleMass/simoptions.maxAllowedStatorModuleMass) + ...
                     1 * design.BaseScore * (10 * design.StatorModuleMass/simoptions.maxAllowedStatorModuleMass)^2;

                score = score + design.maxAllowedStatorModuleMassPenalty;
                
            end
        end
    end
    
    
end



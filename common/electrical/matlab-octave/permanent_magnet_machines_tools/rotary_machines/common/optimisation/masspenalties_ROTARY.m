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
                    10 * design.OptimInfo.BaseScore * (design.RotorMass/simoptions.maxAllowedRotorMass) + ...
                     1 * design.OptimInfo.BaseScore * (10 * design.RotorMass/simoptions.maxAllowedRotorMass)^2;

                score = score + design.maxAllowedRotorMassPenalty;
                
            end
        end
    else
        [score, design, simoptions] = ...
            addpenalty_AM(design, simoptions, 'upper', 'RotorMass', design.OptimInfo.BaseScore, score);
    end

    design.maxAllowedRotorModuleMassPenalty = 0;

    % exceeding max allowed rotor module mass
    if isfield(simoptions, 'maxAllowedRotorModuleMass')
        if ~isempty(simoptions.maxAllowedRotorModuleMass)
            if design.RotorModuleMass > simoptions.maxAllowedRotorModuleMass
                
                design.maxAllowedRotorModuleMassPenalty = ...
                    10 * design.OptimInfo.BaseScore * (design.RotorModuleMass/simoptions.maxAllowedRotorModuleMass) + ...
                     1 * design.OptimInfo.BaseScore * (10 * design.RotorModuleMass/simoptions.maxAllowedRotorModuleMass)^2;

                score = score + design.maxAllowedRotorModuleMassPenalty;
                
            end
        end
    else
        [score, design, simoptions] = ...
            addpenalty_AM(design, simoptions, 'upper', 'RotorModuleMass', design.OptimInfo.BaseScore, score);
    end

    % exceeding max stator mass
    design.maxAllowedStatorMassPenalty = 0;

    if isfield(simoptions, 'maxAllowedStatorMass')
        if ~isempty(simoptions.maxAllowedStatorMass)
            if design.StatorMass > simoptions.maxAllowedStatorMass
                
                design.maxAllowedStatorMassPenalty = ...
                    10 * design.OptimInfo.BaseScore * (design.StatorMass/simoptions.maxAllowedStatorMass) + ...
                     1 * design.OptimInfo.BaseScore * (10 * design.StatorMass/simoptions.maxAllowedStatorMass)^2;

                score = score + design.maxAllowedStatorMassPenalty;
                
            end
        end
    else
        [score, design, simoptions] = ...
            addpenalty_AM(design, simoptions, 'upper', 'StatorMass', design.OptimInfo.BaseScore, score);
    end

    design.maxAllowedStatorModuleMassPenalty = 0;

    % exceeding max stator module mass
    if isfield(simoptions, 'maxAllowedStatorModuleMass')
        if ~isempty(simoptions.maxAllowedStatorModuleMass)
            if design.StatorModuleMass > simoptions.maxAllowedStatorModuleMass
                
                design.maxAllowedStatorModuleMassPenalty = ...
                    10 * design.OptimInfo.BaseScore * (design.StatorModuleMass/simoptions.maxAllowedStatorModuleMass) + ...
                     1 * design.OptimInfo.BaseScore * (10 * design.StatorModuleMass/simoptions.maxAllowedStatorModuleMass)^2;

                score = score + design.maxAllowedStatorModuleMassPenalty;
                
            end
        end
    else
        [score, design, simoptions] = ...
            addpenalty_AM(design, simoptions, 'upper', 'StatorModuleMass', design.OptimInfo.BaseScore, score);
    end
    
    
    % total mass constraint
    [score, design, simoptions] = ...
            addpenalty_AM(design, simoptions, 'upper', 'TotalMass', design.OptimInfo.BaseScore, score);
    
    
end



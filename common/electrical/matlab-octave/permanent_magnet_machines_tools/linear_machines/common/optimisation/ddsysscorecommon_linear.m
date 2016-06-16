function [score, design] = ddsysscorecommon_linear(design, simoptions)
% calculates aspects of the score for a heaving buoy system attached to a
% heaving buoy common to all linear machines.
%
% Syntax:
%
% [score, design] = ddsysscorecommon_linear(design, simoptions)
%
% Input:
% 
% design and simoptions are structures which contain information about the
% machine design and the simulation scoring options. The design structure
% must contain the following fields:
%
%   PowerLoadMean: the mean power exported to the grid by the system
%
%   CostEstimate: an estimate of the cost the machine used in one
%   unit of the system.
%
% However, if other scoring criteria are supplied in the simoptions
% structure, other fields may also be required. The possible scoring
% criteria, and related design structure fields are now summarized. The
% simoptions field, and associated design structure are denoted by 
%
% simfield [designfield]
% 
% If design field is empty (e.g. []) there is no associated additional
% design field.
%
% 1. FarmSize [], specifies that the machines are part of an array to be of
% total generation capacity 'Farmsize'.
%
% 2. MinDesiredMeanPower [], specifies a minimum desired power output per
% system unit. This will be compared to design.GridMeanPower
%
% 3. 
%
%
% Output:
%
% score - the resulting score assigned to the design based on the inputs
% 
% design - machine design matrix returned with a breakdown of any penalties
% applied to the design. The pre-penalty machine score is appended in the
% field 'CostScore'
%
    
% Copyright Richard Crozer, The University of Edinburgh

    % use the amortised cost per kW hour as the base score
    [score, design] = costscore_AM(design, simoptions);
    
     % store the base score in the design structure
    design.BaseScore = score;
    
    % determine any electrical penalties
    [score, design] = electricalpenalties_AM(design, simoptions, score);
    
    % penalty 1: exceeding +/- a certain displacement in the armature

    design.maxAllowedxApenalty = 0;

    if isfield(simoptions, 'maxAllowedxA')
        if ~isempty(simoptions.maxAllowedxA)

            if design.xAmax > simoptions.maxAllowedxA

                design.maxAllowedxApenalty = ...
                    5 * design.BaseScore * (design.xAmax/simoptions.maxAllowedxA);

                score = score + design.maxAllowedxApenalty;

            end
        end
    end

%     design.maxAllowedxTpenalty = 0;
% 
%     if isfield(simoptions, 'maxAllowedxTpenalty')
%         if ~isempty(simoptions.maxAllowedxT)
% 
%             if design.xTmax > simoptions.maxAllowedxT
% 
%                 design.maxAllowedxTpenalty = ...
%                       5 * design.BaseScore * (design.xTmax/simoptions.maxAllowedxT);
% 
%                 score = score + design.maxAllowedxTpenalty;
% 
%             end
%         end
% 
%     end

    % exceeding the mass of half the water displaced by the buoy
    if isfield(simoptions, 'BuoyParameters')
        
        design.buoyMassPenalty = 0;

        maxmass = simoptions.BuoySim.BuoyParameters.mass_external / 2;  
        
        if ~isfield(simoptions, 'NoOfMachines')
            simoptions.NoOfMachines = 1;
        end
        
        totmass = design.massT * simoptions.NoOfMachines;

        if totmass > maxmass
            
            design.buoyMassPenalty = ...
                100*(totmass/maxmass) + (100 * (totmass/maxmass))^2;
            
            score = score + design.buoyMassPenalty;

        end

    end

end
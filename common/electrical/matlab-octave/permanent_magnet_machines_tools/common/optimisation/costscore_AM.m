function [score, design] = costscore_AM(design, simoptions)

    if ~isinf(design.PowerLoadMean)
        
        if all(isfield(simoptions.Evaluation, {'DiscountRate', 'ProjectYears', 'CapacityFactor'}))
            
            % if a specified wave farm size has been given, determine the score
            % based on this
            if isfield(simoptions, 'FarmSize')

                % calculate the required number of devices to make up the
                % desired farm size
                design.RequiredDevices = ceil(simoptions.FarmSize / design.PowerLoadMean);
                % Estimate the total cost of the wave farm
                design.FarmCostEstimate = design.CostEstimate * design.RequiredDevices;
                % calculate the amortized cost in euros per kWhr
                design.CostPerkWhr = amortizedenergycostperkwhr(design.FarmCostEstimate, ...
                                                         design.RequiredDevices * design.PowerLoadMean, ...
                                                         simoptions.Evaluation.DiscountRate, ...
                                                         simoptions.Evaluation.ProjectYears, ...
                                                         simoptions.Evaluation.CapacityFactor);

            else

                % calculate the amortized cost in euros per kWhr
                design.CostPerkWhr = amortizedenergycostperkwhr(design.CostEstimate, ...
                                                         design.PowerLoadMean, ...
                                                         simoptions.Evaluation.DiscountRate, ...
                                                         simoptions.Evaluation.ProjectYears, ...
                                                         simoptions.Evaluation.CapacityFactor);
            end
        
        else
            % calculate the cost per kWhr of output
            design.CostPerkWhr = design.CostEstimate ./ design.PowerLoadMean;
        end
        
    else
        score = Inf;
        return;
    end
    
    % multiply the amortized cost by a scale factor (this is used
    % to change to pence per kWhr from pounds per kWHr for example)
    % the default scale factor is 100
    score = design.CostPerkWhr .* simoptions.Evaluation.CostScaleFactor;
    
end
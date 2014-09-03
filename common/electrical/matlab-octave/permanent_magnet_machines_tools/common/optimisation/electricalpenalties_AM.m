function [score, design, simoptions] = electricalpenalties_AM(design, simoptions, score)
% evaluates penalty functions related to electrical aspects of a generator
%
% Syntax
%
% [score, design, simoptions] = electricalpenalties_AM(design, simoptions)
% [score, design, simoptions] = electricalpenalties_AM(..., score)
%
% Description
%
% electricalpenalties_AM generates penalties using information provided in
% design and simoptions which are structures with fields containing the
% necessary information to evaluate the penalty functions.
% 
% The penalties which are evaluated are determined by the presence or
% absence of fields in the simoptions structure. The information necessary
% to evaluate the penalty must also be present in the design structure. The
% possible penalties and required conbinations of fields are detailed
% below.
%
% 
% Max RMS Current density in wire - Penalty for exceeding a desired maximum
%   current density in the conductor. The field 'maxAllowedJrms' must be
%   present in simoptions. This should contain a scalar value of the
%   maximum desired rms current density in the conductor. If supplied, the
%   field 'JCoilRms' must be present in the design structure. This should
%   contain the rms current in the conductor. The field
%   'maxAllowedJrmspenalty' is added to the design structure containing the
%   associated penalty. This will be zero if the maximum rms current
%   density is not exceeded.
%
% Max Peak Current Density in Wire - Penalty for exceeding a desired peak
%   currrent density in the conductor. The field 'maxAllowedJpeak' must be
%   present in simoptions. This should contain a scalar value of the
%   maximum desired peak current density in the conductor. If supplied, the
%   field 'JCoilPeak' must be present in the design structure. This should
%   contain the peak current density in the conductor. The field
%   'maxAllowedJrmspenalty' is added to the design structure containing the
%   associated penalty. This will be zero if the maximum peak current
%   density is not exceeded.
% 
% Max Peak EMF -
%
% Min RMS EMF - 
%
% Target Mean Power - 
%

    if nargin < 3
        score = 0;
    end

    % exceeding desired maximum rms current desity
    
    % support legacy code with old penalty name
    if isfield(simoptions, 'maxAllowedJrms') ...
            && ~isfield(simoptions, 'max_JCoilRms')
        
        simoptions.max_JCoilRms = simoptions.maxAllowedJrms;
        
        if isfield(simoptions, 'maxAllowedJrmsPenFactor')
            simoptions.max_JCoilRms_penfactor = simoptions.maxAllowedJrmsPenFactor;
        end
    end
    
    % apply the target mean load power penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'JCoilRms', design.OptimInfo.BaseScore, score);

    design.maxAllowedJpeakpenalty = 0;

    % exceeding max allowed peak current density
    
    % support legacy code with old penalty name
    if isfield(simoptions, 'maxAllowedJpeak') ...
            && ~isfield(simoptions, 'max_JCoilPeak')
        
        simoptions.max_JCoilPeak = simoptions.maxAllowedJpeak;
        
        if isfield(simoptions, 'maxAllowedJpeakPenFactor')
            simoptions.max_JCoilPeak_penfactor = simoptions.maxAllowedJpeakPenFactor;
        end
    end
    
    % apply the upper coil peak current density
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'JCoilPeak', design.OptimInfo.BaseScore, score);

    % exceeding max allowed peak phase EMF
    % support legacy code with old penalty name
    if isfield(simoptions, 'maxAllowedEMFpeak') ...
            && ~isfield(simoptions, 'max_EMFPhasePeak')
        
        simoptions.max_EMFPhasePeak = simoptions.maxAllowedEMFpeak;
        
        if isfield(simoptions, 'maxAllowedEMFpeakPenFactor')
            simoptions.max_EMFPhasePeak_penfactor = simoptions.maxAllowedEMFpeakPenFactor;
        end
    end
    
    % apply the target mean load power penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'EMFPhasePeak', design.OptimInfo.BaseScore, score);
    
    % exceeding max allowed rms voltage
    % support legacy code with old penalty name
    if isfield(simoptions, 'maxAllowedRMSEMF') ...
            && ~isfield(simoptions, 'max_EMFPhaseRms')
        
        simoptions.max_EMFPhaseRms = simoptions.maxAllowedRMSEMF;
        
        if isfield(simoptions, 'maxAllowedRMSEMFPenFactor')
            simoptions.max_EMFPhaseRms_penfactor = simoptions.maxAllowedRMSEMFPenFactor;
        end
    end
    
    % apply the target mean load power penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'EMFPhaseRms', design.OptimInfo.BaseScore, score);

    % difference from minimum desired rms emf 
    % support legacy code with old penalty name
    if isfield(simoptions, 'minAllowedRMSEMF') ...
            && ~isfield(simoptions, 'min_EMFPhaseRms')
        
        simoptions.min_EMFPhaseRms = simoptions.minAllowedRMSEMF;
        
        if isfield(simoptions, 'minAllowedRMSEMFPenFactor')
            simoptions.min_EMFPhaseRms_penfactor = simoptions.minAllowedRMSEMFPenFactor;
        end
    end
    
    % apply the target mean load power penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'lower', 'EMFPhaseRms', design.OptimInfo.BaseScore, score);
    
    
%     design.minAllowedEMFrmspenalty = 0;
% 
%     if isfield(simoptions, 'minAllowedRMSEMF')
%         if ~isempty(simoptions.minAllowedRMSEMF)
%             if design.EMFPhaseRms < simoptions.minAllowedRMSEMF
%                 
%                 simoptions = setfieldifabsent(simoptions, 'minAllowedRMSEMFPenFactor', [20, 0]);
%                 
%                 if isscalar(simoptions.minAllowedRMSEMFPenFactor)
%                     simoptions.minAllowedRMSEMFPenFactor = [simoptions.minAllowedRMSEMFPenFactor, 0];
%                 end
%                 
%                 design.minAllowedEMFrmspenalty = ...
%                     simoptions.minAllowedRMSEMFPenFactor(1) * design.OptimInfo.BaseScore * (simoptions.minAllowedRMSEMF/design.EMFPhaseRms) ...
%                       + design.OptimInfo.BaseScore * (simoptions.minAllowedRMSEMFPenFactor(2) * simoptions.minAllowedRMSEMF/design.EMFPhaseRms)^2;
% 
%                 score = score + design.minAllowedEMFrmspenalty;
%                 
%             end
%         end
%     end
    
    % target power is desired
    % support legacy code with old penalty name
    if isfield(simoptions, 'targetMeanPower') ...
            && ~isfield(simoptions, 'target_PowerLoadMean')
        simoptions.target_PowerLoadMean = simoptions.targetMeanPower;
        
        if isfield(simoptions, 'targetMeanPowerPenFactor')
            simoptions.target_PowerLoadMean_penfactor = simoptions.targetMeanPowerPenFactor;
        end
    end
    
    % apply the target mean load power penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'target', 'PowerLoadMean', design.OptimInfo.BaseScore, score);

%     design.targetMeanPowerpenalty = 0;
% 
%     if isfield(simoptions, 'targetMeanPower')
%         if ~isempty(simoptions.targetMeanPower)
%             
%             % we will try to hit target power within a band of 10%
%             
%             minpower = 0.95 * simoptions.targetMeanPower;
%             maxpower = 1.05 * simoptions.targetMeanPower;
%             
%             if design.PowerLoadMean < minpower
%                 
%                 simoptions = setfieldifabsent(simoptions, 'minAllowedMeanPowerPenFactor', [10, 0]);
%                 
%                 if isscalar(simoptions.minAllowedMeanPowerPenFactor)
%                     simoptions.minAllowedMeanPowerPenFactor = [simoptions.minAllowedMeanPowerPenFactor, 0];
%                 end
%                 
%                 design.targetMeanPowerpenalty = ...
%                     simoptions.minAllowedMeanPowerPenFactor(1) * design.OptimInfo.BaseScore * (minpower/design.PowerLoadMean) ...
%                       + design.OptimInfo.BaseScore * (simoptions.minAllowedMeanPowerPenFactor(2) * minpower/design.PowerLoadMean)^2;
% 
%                 score = score + design.targetMeanPowerpenalty;
%                 
%             elseif design.PowerLoadMean > maxpower
%                 
%                 simoptions = setfieldifabsent(simoptions, 'maxAllowedMeanPowerPenFactor', [10, 0]);
%                 
%                 if isscalar(simoptions.maxAllowedMeanPowerPenFactor)
%                     simoptions.maxAllowedMeanPowerPenFactor = [simoptions.maxAllowedMeanPowerPenFactor, 0];
%                 end
%                 
%                 design.targetMeanPowerpenalty = ...
%                     simoptions.maxAllowedMeanPowerPenFactor(1) * design.OptimInfo.BaseScore * (design.PowerLoadMean / maxpower) ...
%                       + design.OptimInfo.BaseScore * (simoptions.maxAllowedMeanPowerPenFactor(2) * design.PowerLoadMean / maxpower)^2;
% 
%                 score = score + design.targetMeanPowerpenalty;
%                 
%             end
%         end
%     end
    
    % minimum desired power
    % support legacy code with old penalty name
    if isfield(simoptions, 'minPowerLoadMean') ...
            && ~isfield(simoptions, 'min_PowerLoadMean')
        
        simoptions.min_PowerLoadMean = simoptions.minPowerLoadMean;
        
        if isfield(simoptions, 'minPowerLoadMeanPenFactor')
            simoptions.min_PowerLoadMean_penfactor = simoptions.minPowerLoadMeanPenFactor;
        end
    end
    
    % apply the penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'lower', 'PowerLoadMean', design.OptimInfo.BaseScore, score);
    
%     if isfield(simoptions, 'minPowerLoadMean')
%         if design.PowerLoadMean < simoptions.minPowerLoadMean
% 
%             simoptions = setfieldifabsent(simoptions, 'minPowerLoadMeanPenFactor', [10, 0]);
% 
%             if isscalar(simoptions.minPowerLoadMeanPenFactor)
%                 simoptions.minPowerLoadMeanPenFactor = [simoptions.minPowerLoadMeanPenFactor, 0];
%             end
% 
%             design.minPowerLoadMeanPenalty = ...
%                 simoptions.minPowerLoadMeanPenFactor(1) * design.OptimInfo.BaseScore * (simoptions.minPowerLoadMean/design.PowerLoadMean) ...
%                 + design.OptimInfo.BaseScore * (simoptions.minPowerLoadMeanPenFactor(2) * simoptions.minPowerLoadMean/design.PowerLoadMean)^2;
% 
%             score = score + design.minPowerLoadMeanPenalty;
%         else
%             design.minPowerLoadMeanPenalty = 0;
%         end
%     end
    
    % support legacy code with old penalty name
    if isfield(simoptions, 'maxPowerLoadMean') ...
            && ~isfield(simoptions, 'max_PowerLoadMean')
        
        simoptions.max_PowerLoadMean = simoptions.maxPowerLoadMean;
        
        if isfield(simoptions, 'maxPowerLoadMeanPenFactor')
            simoptions.max_PowerLoadMean_penfactor = simoptions.maxPowerLoadMeanPenFactor;
        end
    end
    
    % apply the penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'PowerLoadMean', design.OptimInfo.BaseScore, score);
    
    % difference from minimum desired power factor
    % support legacy code with old penalty name
    if isfield(simoptions, 'minAllowedPowerFactor') ...
            && ~isfield(simoptions, 'min_PowerFactor')
        
        simoptions.min_PowerFactor = simoptions.minAllowedPowerFactor;
        
        if isfield(simoptions, 'minAllowedPowerFactorPenFactor')
            simoptions.min_PowerFactor_penfactor = simoptions.minAllowedPowerFactorPenFactor;
        end
    end
    
    % apply the penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'lower', 'PowerFactor', 1, score);
    
    
    % Total harmonic distotion (THD) penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'VoltagePercentTHD', 1, score);
    
    
%     design.minAllowedPowerFactorpenalty = 0;
% 
%     if isfield(simoptions, 'minAllowedPowerFactor')
%         if ~isempty(simoptions.minAllowedPowerFactor)
%             if design.PowerFactor < simoptions.minAllowedPowerFactor
%                 
%                 simoptions = setfieldifabsent(simoptions, 'minAllowedPowerFactorPenFactor', [20, 0]);
%                 
%                 if isscalar(simoptions.minAllowedPowerFactorPenFactor)
%                     simoptions.minAllowedPowerFactorPenFactor = [simoptions.minAllowedPowerFactorPenFactor, 0];
%                 end
%                 
%                 design.minAllowedPowerFactorpenalty = ...
%                     simoptions.minAllowedPowerFactorPenFactor(1) * design.OptimInfo.BaseScore * (simoptions.minAllowedPowerFactor/design.PowerFactor) ...
%                       + design.OptimInfo.BaseScore * (simoptions.minAllowedPowerFactorPenFactor(2) * simoptions.minAllowedPowerFactor/design.PowerFactor)^2;
% 
%                 score = score + design.minAllowedPowerFactorpenalty;
%                 
%             end
%         end
%     end
    
    
end
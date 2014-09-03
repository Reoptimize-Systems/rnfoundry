function [score, design, simoptions] = structuralpenalties_AM(design, simoptions, score)
% calculates structural penalties for a design
%
% Syntax
%
% [score, design, simoptions] = structuralpenalties_AM(design, simoptions)
% [score, design, simoptions] = structuralpenalties_AM(design, simoptions, score)
%
% Description
%
% structuralpenalties_AM generates a penalty related to the maximum
% deflection in a structure. The input design and simoptions structures
% contain fields which determine the penalty. 
%
% simoptions can contain the follwoing fields:
%   maxAllowedDeflection - scalar value the maximum permitted deflection,
%     above which a penalty will be applied
%
%   maxAllowedDeflectionFactor - scalar value between 0 and 1 determining
%     how much deflection is permitted in the design being evaluated. The
%     allowed deflection is determined by multiplying this factor by the
%     value supplied in a field 'g' in the design structure, i.e.
%
%     Max Allowed Deflection = maxAllowedDeflectionFactor * design.g
%
%     The value of the maximum allowed deflection is compared with the
%     value of the actual deflection which must be provided in another
%     field in the design structure, 'MaxDeflection'. 
%
%   maxAllowedDeflectionPenFactor - either a scalar value, or a 2 element
%    vector which determines the size of the penalty added if the maximum
%    deflection exceeds the maximum allowed deflection. These represent
%    penalty factors, i.e. [p1, p2]. If a scalar value is supplied this is
%    the p1 value and p2 is set to zero. The resulting penalty is then
%    calculated using the formula
%
%    penalty = (base score * p1 * Max Deflection / Max Allowed Deflection)
%               + (base score * (p2 * Max Deflection / Max Allowed Deflection)^2)
%
%    The base score must be provided in the field 'BaseScore' in the design
%    structure. The value of the penalty is also added to the design
%    strucuture in the field 'maxAllowedDeflectionPenalty'.
%
% An optional third argument is an existing score to which the penalty must
% be added. If not supplied this is set to zer, and the returned value is
% only the penalty as calculated using the supplied parameters. 
% 
% If both maxAllowedDeflection and maxAllowedDeflectionFactor are not
% present in the simoptions structure, no penalty is applied and nothing is
% added to the design structure. In this case the score added is zero.
%
%

% Copyright 2013 Richard Crozier (richard.crozier@yahoo.co.uk)

    if nargin < 3
        score = 0;
    end
    
    % exceeding max allowed rotor mass
    design.maxAllowedDeflectionPenalty = 0;

    if isfield(simoptions, 'maxAllowedDeflection') ...
        && isfield(simoptions, 'maxAllowedDeflectionFactor') 
        
        if ~isempty(simoptions.maxAllowedDeflectionFactor)
            
            % calculate the max allowed deflection from 
            maxAllowedDeflection = simoptions.maxAllowedDeflectionFactor * design.g;
            
            % set an infinite max allowed absolute deflection if it has not
            % already been set
            if isempty(simoptions.maxAllowedDeflection)
                simoptions.maxAllowedDeflection = inf;
            end
            
            % compare the calculaed deflection from the deflection factor
            % to the absolute max deflection value, and replace if it is
            % smaller
            if maxAllowedDeflection < simoptions.maxAllowedDeflection
                simoptions.maxAllowedDeflection = maxAllowedDeflection;
            end
            
            simoptions = setfieldifabsent(simoptions, 'maxAllowedDeflectionPenFactor', [10, 10]);
            
            if isscalar(simoptions.maxAllowedDeflectionPenFactor)
                simoptions.maxAllowedDeflectionPenFactor = [simoptions.maxAllowedDeflectionPenFactor, 0];
            end
                
            if design.MaxDeflection > simoptions.maxAllowedDeflection
                
                design.maxAllowedDeflectionPenalty = ...
                    simoptions.maxAllowedDeflectionPenFactor(1) * design.OptimInfo.BaseScore * (design.MaxDeflection/simoptions.maxAllowedDeflection) + ...
                     1 * design.OptimInfo.BaseScore * (simoptions.maxAllowedDeflectionPenFactor(2) * design.MaxDeflection/simoptions.maxAllowedDeflection)^2;

                score = score + design.maxAllowedDeflectionPenalty;
                
            end
        end
    end

end
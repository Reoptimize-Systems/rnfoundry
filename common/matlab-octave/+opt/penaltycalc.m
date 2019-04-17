function penalties = penaltycalc (variables, penaltyspec, varargin)
% creates penalties based on a base score and constraints
%
% Syntax
%
% penalties = penaltycalc (variables, penaltyspec, options.BaseScore)
%
% Description
%
% penaltycalc calcuates penalties based on the values of structure fields
% or object properties. Optionally these can be related to a 'base' score
% such that improvements in the base score also reduce the penlties.
% Penalties can be based upper or lower bounds on the variable or a target
% value. The upper, lower and target values can optionally themselves be
% fields or properties specified in the variables input.
%
%
% Input
%
%   variables - either a normal matlab structure or an object (class
%     instance) with properties.
% 
%   penaltyspec - structure specifying a set of penalties to be applied 
%     based on the value of the fields or properties of 'variables'. Each
%     field in the penaltyspec structure is itself a structure which
%     defines a single penalty. Each of these nested penalty structures can
%     contain the following fields:
%
%     FieldName : optional field which contains the name of the field or
%       property in 'variables' which is to be tested against the penalty
%       specification. If not present, the field name of the penalty
%       specification in the penaltyspec structure is used instead.
% 
%     Type : string determining the type of penalty to be applied. The
%     penalty is applied based on the value in the 'Limit' field (see
%     below). It can be one of the following values:
%   
%       'upper'  In this case the test is for an upper limit on a quantity.
%                If the value of the quantity to be tested exceeds this
%                value, a penalty will be applied. The value of the limit
%                above which the penalty is applied is in the 'Limit' field
%                of the penaltyspec structure.
%
%                The penalty that is applied is a function with
%                coefficients specified in the PenaltyFactors field (see
%                below). The following function is applied when the value
%                of the quantity exceeds the specified limit:
%
%                penalty = PenaltyFactors(1) ...
%                              * BaseScore 
%                              * (quantity_to_test / upper_limit) ...
%                          + BaseScore 
%                              * (PenaltyFactors(2) * quantity_to_test / upper_limit)^2
%
%                The contents of PenaltyFactors can either be a scalar or a
%                two element vector. If a scalar, the second coefficient is
%                automatically set to zero.
%
%       'lower'  In this case the test is for a lower limit on a quantity. 
%                If the value of the quantity to be tested falls below this
%                value, a penalty will be applied. The value of the limit
%                below which the penalty is applied is in the 'Limit' field
%                of the penaltyspec structure.
%
%                The penalty that is applied is a function with
%                coefficients specified in the PenaltyFactors field (see
%                below). The following function is applied when the value
%                of the quantity exceeds the specified limit:
%
%                penalty = PenaltyFactors(1) ...
%                              * BaseScore ...
%                              * (1 + abs((lower_limit - quantity_to_test)/lower_limit)) ...
%          	               + BaseScore ...
%                              * (PenaltyFactors(2) * (1 + abs((lower_limit - quantity_to_test)/lower_limit)))^2;
%
%                The contents of PenaltyFactors can either be a scalar or a
%                two element vector. If a scalar, the second coefficient is
%                automatically set to zero.
%
%       'target' In this case the test is for a target value for the
%                quantity. This applies an upper and lower limit penalty to
%                the quantity with some tolerance around the target value.
%                
%                By default, a tolerance of +/-10% of the target value is
%                used (i.e. 20% variation). An optional field 'TargetTol'
%                can be supplied which should contain a scalar value. This
%                value is a factor which determines the tolerance in terms
%                of the value of the quantity, e.g
%                   
%                penaltyspec.ExampleFieldName.Type = 'target';
%                penaltyspec.ExampleFieldName.Limit = 100;
%                penaltyspec.ExampleFieldName.TargetTol = 0.05;
%
%                sets a target value of 100 and allows a +/-5% tolerance on
%                the value found in ExampleFieldName in the test structure,
%                i.e. +/-5 in this case.
%
%     Limit : field specifying the limit on the quantity to which a penalty
%       is to be applied. This should be either a scalar value, or a
%       character vector containing the name of another field or property
%       in the variables structure. In the second case, the value in the
%       given field or property is used as the limit.
%
%     PenaltyFactors : field containing the coefficients of the 
%       function to be used to calculate the penalty applied when a
%       specified limit is breached. The possible contents of
%       PenaltyFactors depends on the type of limit to be applied. See the
%       description of the 'Type' field above for more information on
%       appropriate values for the PenaltyFactors field for each penalty
%       type.
%
%     BaseScore : scalar value containing a base score to be used when
%       calculating the penalty, see individual penalty type description in
%       the 'Type' field (see above) to see how this value is used for each
%       penalty type. This value overrides, for only this variable, any
%       global base score for all variables specified using the 'BaseScore'
%       option described below (including the default global base score of
%       1 if this option is not used).
%
% Additional options can be supplied as parameter-value pairs. The
% available options are:
%
%   'BaseScore' - optional base score which will be used as the basis for
%     all applied the penalties, unless another alternative base score is
%     provided for a specific penalty. A function either linear or
%     quadratic will be applied to the base score to generate the penalty
%     value. Default is 1 if not supplied. 
%
% Examples
%
% vars.bananas = 10;
% penaltyspec.bananas.Type = '';
% penaltyspec.bananas.Limit = 8;
% penaltyspec.bananas.PenaltyFactors = 5;
% penalties = penaltycalc (variables, penaltyspec);
%
    
    options.BaseScore = 1;
    
    options = parse_pv_pairs (options, varargin);
    
    if ~(isstruct (variables) || isobject (variables))
        error ('variables must be a structure or an object (class instance)');
    end
    
    pennames = fieldnames (penaltyspec);
    
    penalties.total = 0;
    
    for ind = 1:numel (pennames)
        
        this_penalty = penaltyspec.(pennames{ind});
        
        type = this_penalty.Type;
        
        if ~isfield (this_penalty, 'FieldName')
            this_penalty.FieldName = pennames{ind};
        end
        
        if ~isfield (this_penalty, 'BaseScore')
            this_penalty.BaseScore = options.BaseScore;
        end
        
        if isempty (this_penalty.FieldName)
            
            quantityfname = pennames{ind};
            
            if isstruct (variables)
                
                if ~isfield (variables, quantityfname)
                    error ('FieldName for penalty %s is empty, so tried using penalty fieldname, but this is not present in the variables structure', pennames{ind});
                end
                
            elseif isobject (variables)
                
                if isprop (variables, quantityfname)
                    error ('FieldName for penalty %s is empty, so tried using penalty fieldname, but this is not a property of the supplied object in variables', pennames{ind});
                end
                
            end
        else
            quantityfname = this_penalty.FieldName;
        end
        
        varval = variables.(quantityfname);
        
        penalty_factor = this_penalty.PenaltyFactors;
        
        if ischar (this_penalty.Limit)
            limit = variables.(this_penalty.Limit);
        else
            limit = this_penalty.Limit;
        end

        
        switch type

            case 'upper'

                penfieldname = ['max_', quantityfname, '_penalty'];

                % exceeding max allowed value of quantity
                penalties.(penfieldname) = 0;

                if ~isempty(limit)

                    if varval > limit

                        if isempty (penalty_factor)
                            penalty_factor = [1, 0];
                        end

                        if isscalar(penalty_factor)
                            penalty_factor = [penalty_factor, 0];
                        end

                        penalties.(penfieldname) = apply_upper (penalty_factor, this_penalty.BaseScore, varval, limit);
                        
                        penalties.total = penalties.total + penalties.(penfieldname);

                    end

                end


            case 'lower'

                % penalty related to a quantity being below some desired value

                penfieldname = ['min_', quantityfname, '_penalty'];

                penalties.(penfieldname) = 0;

                if ~isempty(limit)

                    if varval < limit

                        if isempty (penalty_factor)
                            penalty_factor = [1, 0];
                        end

                        if isscalar(penalty_factor)
                            penalty_factor = [penalty_factor, 0];
                        end

                        penalties.(penfieldname) = apply_lower (penalty_factor, this_penalty.BaseScore, varval, limit);

                        penalties.total = penalties.total + penalties.(penfieldname);

                    end

                end


            case 'target'

                % hit a target value, with a given amount of allowed variation

                if ~isfield (this_penalty, 'TargetTol')
                    this_penalty.TargetTol = 0.1;
                end
                
                % allowed variation supplied
                if this_penalty.TargetTol <= 0
                    error('TargetTol for %s must be greater than 0', quantityfname);
                end
                
                if numel(limit) ~= 1
                    error('Target penalty limit for %s was the wrong size (%d)', quantityfname, numel(limit))
                end

                if isempty (penalty_factor)
                    penalty_factor = [1, 0; 1, 0];
                end
                
                if size (penalty_factor, 1) == 1
                    % use same penalty factor for upper and lower
                    penalty_factor = [penalty_factor; penalty_factor];
                elseif size (penalty_factor, 1) > 2
                    error ('Target penalty_factor for %s had more than 2 rows', pennames{ind});
                end
                
                if size (penalty_factor, 2) == 1
                    penalty_factor = [penalty_factor, [0; 0]];
                elseif size (penalty_factor, 2) > 2
                    error ('Target penalty_factor for %s had more than 2 columns', pennames{ind});
                end

                maxlimitfieldname = ['target_max_', quantityfname];

                minlimitfieldname = ['target_min_', quantityfname];

                limit_upper = limit(1) + (this_penalty.TargetTol * abs(limit(1)))/2;

                limit_lower = limit(1) - (this_penalty.TargetTol * abs(limit(1)))/2;

                % apply a penalty for exceeding an upper bound on the
                % target
                if varval > limit_upper
                    
                    penalties.(maxlimitfieldname) = apply_upper (penalty_factor(1,:), this_penalty.BaseScore, varval, limit_upper);
                    
                    penalties.total = penalties.total + penalties.(maxlimitfieldname);
                    
                end

                % apply a penalty for falling below a lower bound on the
                % target
                if varval < limit_lower
                    
                    penalties.(minlimitfieldname) = apply_lower (penalty_factor(2,:), this_penalty.BaseScore, varval, limit_lower);
                    
                    penalties.total = penalties.total + penalties.(minlimitfieldname);
                    
                end

            otherwise
                
                error ('Penalty %s contained an invalid penalty type: ''%s''', ...
                       pennames{ind}, type);

        end
        
    end

end


function pen = apply_upper (penalty_factor, base_score, var, limit)

    pen = penalty_factor(1) * base_score * (var / limit)  ...
                             + base_score * (penalty_factor(2) * var / limit)^2;
                         
end

function pen = apply_lower (penalty_factor, base_score, var, limit)

    quantdiff = limit - var;
    
    pen = penalty_factor(1) * base_score * (1 + abs(quantdiff/limit)) ...
          	+ base_score * (penalty_factor(2) * (1 + abs(quantdiff/limit)))^2;
                         
end



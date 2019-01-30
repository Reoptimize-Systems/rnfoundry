function penalties = penaltycalc (variables, penaltyspec, base_score)
% function used to create penalties based on a base score and constraints
%
% Syntax
%
% 
%
% Description
%
% penaltycalc searches a structure of object for the field or proerty name
% specified in the quantityfname variable. Depending on what type of
% penalty is desired (specified in 'type') the simoptions sturture is then
% searched for a two matching fields such as max_<quantity name>,
% min_<quantity name> or target_<quantity name> and max_<quantity
% name>_penfactor etc. which determines what penlty function is applied
% (see documentation for type below for more information). An an example,
% for a maximum limit on a field called bannanas in the design structure it
% works in the following way. If the simoptions field is a single number a
% linear penalty funciton will be appplied calculated like so:
%
%   penalty = penaltyspec.bananas.penalty_factor * base_score * (variables.(penaltyspec.bananas.FieldName) / penaltyspec.bananas.limit) 
% 
% if simoptions is a two element vector a quadratic funtion is applied like
% so:
%
%   penalty = simoptions.max_bananas_penfactor(1) * base_score * (design.bananas / simoptions.max_bananas) ...
%              + base_score * (simoptions.max_bananas_penfactor(2) * design.bananas / simoptions.max_bananas)^2
%
% This would be done like so:
%
% vars.bananas = 10
% penaltyspec.bananas.Limit = 8
% penaltyspec.bananas.PenaltyFactors = 5
% penalties = penaltycalc (variables, penaltyspec)
%
% Input
%
%   variables - either a normal matlab structure or an object (class
%     instance)
% 
%   penaltyspec - structure defining a set of penalties to be applied based
%     on the value of the fields or properties of variables.
% 
%     type : string determining the type of penalty to be applied, one of the
%       following values:
%   
%       'upper'  In this case the test is for an upper limit on a value. If
%                the value of the quantity to be tested exceeds this value, a
%                penalty will be applied. In this case it is expected that
%                there will be a field in the simoptions structure named
%                max_<quantity name>, where <quantity name> is the string
%                provided in the quantityfname input. If the field
%                max_<quantity name>_penfactor is also present in the
%                simoptions structure it will be used to determine the
%                severit of the penalty
%
%       'lower'  In this case the test is for a lower limit on a value. If
%                the value of the quantity to be tested falls below this
%                value, a penalty will be applied.
%
%       'target' In this case the test is for a target value forthe quantity.
%                This applies a max and min penalty to the quantity with some
%                tolerance around the target value. If
%                simoptions.target_<quantity name> is a single value a
%                default tolerance of +/-10% of the target value is used
%                (i.e. 20% variation). If simoptions.target_<quantity name>
%                is a two-element vector, the second value is a factor which
%                determines the tolerance in terms of the value of the
%                quantity, e.g 
%
%                simoptions.target_<quantity name> = [100, 0.05] 
%
%                sets a target value of 100 and allows a +/-5% tolerance on
%                the value, i.e. +/-5 in this case.
%     
% 
%   base_score - optional base score which will be used as the basis for the
%     penalty. A function either linear or quadratic will be applied to the
%     base score to generate the penalty value. Default is 1 if not
%     supplied.
%
% 
    
    if nargin < 3
        base_score = 1;
    end
    
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
        
        limit = this_penalty.Limit;

        
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

                        penalties.(penfieldname) = apply_upper (penalty_factor, base_score, varval, limit);
                        
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

                        penalties.(penfieldname) = apply_lower (penalty_factor, base_score, varval, limit);

                        penalties.total = penalties.total + penalties.(penfieldname);

                    end

                end


            case 'target'

                % hit a target value, with a given amount of allowed variation

                if numel(limit) == 1
                    % allow a variation of 10 % around target
                    limit(2) = 0.1;
                elseif numel(limit) == 2
                    % allowed variation supplied
                    if limit(2) <= 0
                        error('Allowed target variation must be greater than 0');
                    end
                else
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

                limit_upper = limit(1) + (limit(2) * abs(limit(1)))/2;

                limit_lower = limit(1) - (limit(2) * abs(limit(1)))/2;

                % apply a penalty for exceeding an upper bound on the
                % target
                if varval > limit_upper
                    
                    penalties.(maxlimitfieldname) = apply_upper (penalty_factor(1,:), base_score, varval, limit_upper);
                    
                    penalties.total = penalties.total + penalties.(maxlimitfieldname);
                    
                end

                % apply a penalty for falling below a lower bound on the
                % target
                if varval < limit_lower
                    
                    penalties.(minlimitfieldname) = apply_lower (penalty_factor(2,:), base_score, varval, limit_lower);
                    
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



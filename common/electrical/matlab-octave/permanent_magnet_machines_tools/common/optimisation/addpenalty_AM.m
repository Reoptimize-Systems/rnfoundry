function [score, design, simoptions] = addpenalty_AM(design, simoptions, type, quantityfname, base_score, score)
% function used to create penalties and add them to a design structure
%
% Syntax
%
% [score, design, simoptions] = addpenalty_AM(design, simoptions, type, quantityfname, base_score, score)
%
% Description
%
% addpenalty_AM searches the design structure for the field name specified
% in the quantityfname variable. Depending on what type of penalty is
% desired (specified in 'type') the simoptions sturture is then searched
% for a two matching fields such as max_<quantity name>, min_<quantity
% name> or target_<quantity name> and max_<quantity name>_penfactor etc.
% which determines what penlty function is applied (see documentation for
% type below for more information). An an example, for a maximum limit on a
% field called bannanas in the design structure it works in the following
% way. If the simoptions field is a single number a linear penalty funciton
% will be appplied calculated like so:
%
%   penalty = simoptions.max_bananas_penfactor * base_score * (design.bananas / simoptions.max_bananas) 
% 
% if simoptions is a two element vector a quadratic funtion is applied like
% so:
%
%   penalty = simoptions.max_bananas_penfactor(1) * base_score * (design.bananas / simoptions.max_bananas) ...
%              + base_score * (simoptions.max_bananas_penfactor(2) * design.bananas / simoptions.max_bananas)^2
%
% This would be done like so:
%
% score = 0
% base_score = 1
% design.bananas = 10
% simoptions.max_bananas = 8
% simoptions.max_bananas_penfactor = 5
% [score, design, simoptions] = addpenalty_AM(design, simoptions, 'upper', 'bananas', base_score, score)
%
% Input
%
%   design, 
% 
%   simoptions - 
% 
%   type - string determining the type of penalty to be applied, one of the
%     following values:
%   
%     'upper'  In this case the test is for an upper limit on a value. If
%              the value of the quantity to be tested exceeds this value, a
%              penalty will be applied. In this case it is expected that
%              there will be a field in the simoptions structure named
%              max_<quantity name>, where <quantity name> is the string
%              provided in the quantityfname input. If the field
%              max_<quantity name>_penfactor is also present in the
%              simoptions structure it will be used to determine the
%              severit of the penalty
%
%     'lower'  In this case the test is for a lower limit on a value. If
%              the value of the quantity to be tested falls below this
%              value, a penalty will be applied.
%
%     'target' In this case the test is for a target value forthe quantity.
%              This applies a max and min penalty to the quantity with some
%              tolerance around the target value. If
%              simoptions.target_<quantity name> is a single value a
%              default tolerance of +/-10% of the target value is used
%              (i.e. 20% variation). If simoptions.target_<quantity name>
%              is a two-element vector, the second value is a factor which
%              determines the tolerance in terms of the value of the
%              quantity, e.g 
%
%              simoptions.target_<quantity name> = [100, 0.05] 
%
%              sets a target value of 100 and allows a +/-5% tolerance on
%              the value, i.e. +/-5 in this case.
%     
% 
%   quantityfname - name of the field in the design structure containing
%     the quantity which will be compared to the limit to generate the
%     penalty.
% 
%   base_score - base score which will be used as the basis for the
%     penalty. A function either linear or quadratic will be applied to the
%     base score to generate the penalty value. 
%
%   score - the input score plus the newly generated penalty
%
% 

    if nargin < 6
        score = 0;
    end
    
    if nargin < 5
        base_score = 1;
    end

    switch type
        
        case 'upper'
            
            limitfieldname = ['max_', quantityfname];
            penfactfieldname = ['max_', quantityfname, '_penfactor'];
            penfieldname = ['max_', quantityfname, '_penalty'];
            
            % exceeding max allowed value of quantity
            design.OptimInfo.(penfieldname) = 0;

            if isfield(simoptions, limitfieldname)
                
                if ~isempty(simoptions.(limitfieldname))
                    
                    if design.(quantityfname) > simoptions.(limitfieldname)

                        simoptions = setfieldifabsent(simoptions, penfactfieldname, [1, 0]);

                        if isscalar(simoptions.(penfactfieldname))
                            simoptions.(penfactfieldname) = [simoptions.(penfactfieldname), 0];
                        end

                        design.OptimInfo.(penfieldname) = ...
                            simoptions.(penfactfieldname)(1) * base_score * (design.(quantityfname) / simoptions.(limitfieldname))  ...
                             + base_score * (simoptions.(penfactfieldname)(2) * design.(quantityfname) / simoptions.(limitfieldname))^2;

                        score = score + design.OptimInfo.(penfieldname);

                    end
                    
                end
                
            end

        case 'lower'
            
            % penalty related to a quantity being below some desired value
            
            limitfieldname = ['min_', quantityfname];
            penfactfieldname = ['min_', quantityfname, '_penfactor'];
            penfieldname = ['min_', quantityfname, '_penalty'];
    
            design.OptimInfo.(penfieldname) = 0;

            if isfield(simoptions, limitfieldname)
                
                if ~isempty(simoptions.(limitfieldname))
                    
                    if design.(quantityfname) < simoptions.(limitfieldname)

                        simoptions = setfieldifabsent(simoptions, penfactfieldname, [1, 0]);

                        if isscalar(simoptions.(penfactfieldname))
                            simoptions.(penfactfieldname) = [simoptions.(penfactfieldname), 0];
                        end

                        quantdiff = simoptions.(limitfieldname) - design.(quantityfname);
                        
                        design.OptimInfo.(penfieldname) = ...
                            simoptions.(penfactfieldname)(1) * base_score * (1 + abs(quantdiff/simoptions.(limitfieldname))) ...
                              + base_score * (simoptions.(penfactfieldname)(2) * (1 + abs(quantdiff/simoptions.(limitfieldname))))^2;

                        score = score + design.OptimInfo.(penfieldname);

                    end
                    
                end
                
            end

        case 'target'
            
            % hit a target value, with a given amount of allowed variation
            targetfieldname = ['target_', quantityfname];
            targetpenfactfieldname = ['target_', quantityfname, '_penfactor'];
            
            if isfield(simoptions, targetfieldname) && ~isempty(simoptions.(targetfieldname))
                
                if numel(simoptions.(targetfieldname)) == 1
                    % allow a variation of 10 % around target
                    simoptions.(targetfieldname)(2) = 0.1;
                elseif numel(simoptions.(targetfieldname)) == 2
                    % allowed variation supplied
                    if simoptions.(targetfieldname)(2) <= 0;
                        error('Allowed target variation must be greater than 0');
                    end
                else
                    error('simoptions.%s was the wrong size (%d)', targetfieldname, numel(simoptions.(targetfieldname)))
                end
                
                if ~isfield(simoptions, targetpenfactfieldname)
                    simoptions.(targetpenfactfieldname) = [1,0];
                end

                maxlimitfieldname = ['max_', quantityfname];
                maxpenfactfieldname = ['max_', quantityfname, '_penfactor'];

                minlimitfieldname = ['min_', quantityfname];
                minpenfactfieldname = ['min_', quantityfname, '_penfactor'];

                simoptions.(maxlimitfieldname) = simoptions.(targetfieldname)(1) ...
                            + (simoptions.(targetfieldname)(2) * abs(simoptions.(targetfieldname)(1)))/2;

                simoptions.(minlimitfieldname) = simoptions.(targetfieldname)(1) ...
                            - (simoptions.(targetfieldname)(2) * abs(simoptions.(targetfieldname)(1)))/2;
                        
                simoptions.(maxpenfactfieldname) = simoptions.(targetpenfactfieldname);
                simoptions.(minpenfactfieldname) = simoptions.(targetpenfactfieldname);
                
                % apply a penalty for exceeding an upper bound on the
                % target
                [score, design, simoptions] = addpenalty_AM(design, simoptions, 'upper', quantityfname, base_score, score);
                
                % apply a penalty for falling below a lower bound on the
                % target
                [score, design, simoptions] = addpenalty_AM(design, simoptions, 'lower', quantityfname, base_score, score);
            
            end
            
        otherwise
            
    end

end
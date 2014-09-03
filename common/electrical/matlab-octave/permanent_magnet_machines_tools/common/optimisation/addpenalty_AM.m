function [score, design, simoptions] = addpenalty_AM(design, simoptions, type, quantityfname, base_score, score)
% function used to create penalties and add them to a design structure
%
% Syntax
%
% [score, design, simoptions] = addpenalty_AM(design, simoptions, type, quantityfname, base_score, score)
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
%     'target' In this case the test is for a lower limit on a value. If
%              the value of the quantity to be tested falls below this
%              value, a penalty will be applied.
%
%     
% 
%   quantityfname - name of the field in the design structure containing
%     the quantity which will be compared to the limit to generate the
%     penalty.
% 
%   score
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
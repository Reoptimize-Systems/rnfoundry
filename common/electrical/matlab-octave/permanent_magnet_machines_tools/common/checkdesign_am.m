function design = checkdesign_am(design, valfields, ratiofields, basefield, varargin)
% checkdesign_am checks a design structure for consistancy given a set
% dimension names and associated ratio names, and can also check that
% ratios are within a given range
%
% Syntax
% 
% design = checkdesign_am(design, valfields, ratiofields, basefield)
% design = checkdesign_am(..., 'Parameter', 'Value')
%
% Input
%
% design is a structure containing fields which represent either absolute
%   values or ratios. It must contain the field name passed in in
%   'basefield'.
%
% valfields is a cell array of strings containing the names of the minimum
%   set of absolute value fields in the sturcture which must be present in
%   order for it to be considered complete.
%
% ratiofields is a cell array of strings containing the names of the
%   minimum set of ratios, the numerators and denominators of which are
%   made up of combinations of the minimum set of value fields in valfields
%   and the base value field in basefield.
%
% basefield is a string containing the name of a field in the structure
%   to which the fields in ratiofields are ultimately based. i.e. all of
%   the dimension fields should be able to be calculated provided this base
%   value is available. 
%
% Description
% 
% design = checkdesign_am(design, valfields, ratiofields, basefield) looks
% in the structure design for field names in the cell arrays valfields,
% ratiofields and the string basefield. If all of the fields are found, it
% checks that the ratio fields are consistant with the corresponding value
% fields. If they are not consistant checkdesign_am will by default throw
% an error. Alternatively the precendence of either the value fields or the
% ratio fields can be chosen through the parameter-value pair 'Precedence',
% where the value must be either the string 'ratios' or 'values'. In this
% case either the values or ratios are overwritten to be consistent with
% the set specified by this parameter.
%
% If only one complete set of either the value fields or ratio fields are
% found, the incomplete set is overwritten and completed to be consistent
% with the complete set.
%
% After completing the field sets, checkdesign_am then checks the ratio
% fields for appropriate values. By default is is assumed that all ratios
% must be positive. This can be overridden through the 'Positive' p-v pair,
% where the value can either be a single boolean, where true indicates all
% ratios must be positive, or false indicates no ratios must be positive,
% or it can be a matrix of booleans, each one corresponding to a ratio,
% determining whether that ratio must be positive.
%
% Similarly a range bounding each ratio can be specified using the
% 'RatioLimits' p-v pair where the value if an (n x 2) matrix of values,
% each row denoting the minimum and maximum values allowed for the
% corresponding ratio.
%
% 


    Inputs.Positive = true;
    Inputs.RatioLimits = [];
    Inputs.ZeroToOne = [];
    Inputs.ValuesTests = {}; % some day we could use this to pass in function handles for value tests
    Inputs.Precedence = '';
    Inputs.RatioSplitString = 'V';
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    if isempty(Inputs.ZeroToOne)
        Inputs.ZeroToOne = repmat(false, 1, numel(ratiofields));
    end
    
    if numel(Inputs.Positive) == 1
        Inputs.Positive = repmat(Inputs.Positive, 1, numel(ratiofields));
    elseif numel(Inputs.Positive) ~= numel(ratiofields)
        error('CHECKDESIGN_AM:badratiolimits', ...
             ['If supplying more than one Positive specification value, ', ...
              'you must supply the same number of values as there are ratiofields.']);
    end
    
    if ~isempty(Inputs.RatioLimits) && ~all(size(Inputs.RatioLimits) == [numel(ratiofields), 2])
        error('CHECKDESIGN_AM:badratiolimits', ...
            ['If limits for the ratios are supplied they must be a two ',...
             'column matrix with the same number of rows as the number ',...
             'of supplied ratios.'])
    end
    
    % check to see if we have a full set of either the required dimension
    % fields or ratio fields
    hasfields = checkforfieldsets(design, valfields, ratiofields);
    
    if all(hasfields)
        
        % get the numerators and denominators from the ratio fields
        [numerstr, denomstr] = splitratiostr(ratiofields, Inputs.RatioSplitString);
        
        % check they all exist
        if ~all(checkforfieldsets(design, unique(numerstr), unique(denomstr)))
           error('CHECKDESIGN_AM:numersdenomsanddimsnotmatch', ...
               ['The numerators and denominators from the split ratio fields ',...
                'could not all be found in the design structure, although ',...
                'both minimum sets of value and ratio fields specified ',...
                'were present']) 
        end
        
        % check if the dimensions and ratios present in the structure are
        % consistant
        
        temp1 = structvals2structratios(design, numerstr, denomstr, Inputs.RatioSplitString);
        
        temp2 = structratios2structvals(design, ratiofields, base, Inputs.RatioSplitString);
        
        ratiosandvalsmatch = true;
        
        for i = 1:numel(ratiofields) 
            
            if temp1.(ratiofields{i}) ~= design.(ratiofields{i})
                
                ratiosandvalsmatch = false;
                
                break;
                
            end
            
        end
        
        if ratiosandvalsmatch || temp2.(basefield) ~= design.(basefield) || temp1.(basefield) ~= design.(basefield)
            
            for i = 1:numel(fields1)

                if temp2.(fields1{i}) ~= design.(fields1{i}) || temp2.(fields2{i}) ~= design.(fields2{i})

                    ratiosandvalsmatch = false;
                    
                    break;

                end

            end
            
        else
            
            ratiosandvalsmatch = false;
            
        end
        
        % clear the temporary structures used for the check
        clear('temp1');
        clear('temp2');
        
        if ~isempty(Inputs.Precedence)
            
            switch Inputs.Precedence
                
                case 'ratios'
                    
                    warning('CHECKDESIGN_AM:ratioprecendence', ...
                            ['The ratio fields do not match the value fields. You\n', ...
                             'have specified that ratios take precedence over the \n',...
                             'dimensions, so the dimensions will be updated to reflect \n',...
                             'the ratios']);
                         
                     design = structratios2structvals(design, ratiofields, base, Inputs.RatioSplitString);
                    
                case 'values'
                    
                    warning('CHECKDESIGN_AM:dimprecendence', ...
                            ['The value fields do not match the ratio fields. You\n', ...
                             'have specified that dimensions take precedence over the \n',...
                             'ratios, so the ratios will be updated to reflect \n',...
                             'the dimensions.']);
                    
                otherwise
                    
                    error('CHECKDESIGN_AM:badprecendencestring', ...
                        'Unknown value\\ratio precendence specification string.')
                    
            end
            
        else
            error('CHECKDESIGN_AM:noprecendence', ...
                 ['The ratio fields do not match the dimension fields and you ',...
                  'have not specified which should take precendence in this case.']);
        end
        
    elseif hasfields(1)
        
        % get the numerators and denominators from the ratio fields
        [numerstr, denomstr] = splitratiostr(ratiofields, Inputs.RatioSplitString);
        
        % check they all exist
        if ~all(checkforfieldsets(design, unique(numerstr), unique(denomstr)))
           error('CHECKDESIGN_AM:numersdenomsanddimsnotmatch', ...
               ['The numerators and denominators from the split ratio fields ',...
                'could not all be found in the design structure, although ',...
                'the minimum set of specified dimension fields are present.',...
                'were present']) 
        end
        
        
        design = structvals2structratios(design, numerstr, denomstr, Inputs.RatioSplitString);
        
        
    elseif hasfields(2)
        
        design = structratios2structvals(design, ratiofields, basefield, Inputs.RatioSplitString);
        
    else
        error('CHECKDESIGN_AM:notenoughfields', ...
            'There is neither a full set of dimension fields of ratio fields.')
    end
    
    % Now check the ratios
    ratioresult = checkratios(design, ratiofields, Inputs.Positive, Inputs.ZeroToOne, Inputs.RatioLimits);
    
    if anyalldims(ratioresult) 
        error('CHECKDESIGN_AM:badratios', ...
            'Ratios did not pass specified tests');
    end
        
end


function ratioresult = checkratios(design, ratiofields, positive, zerotoone, ratiolimits)
% checks that each of a set of ratios in a structure meets various
% specified criteria

    ratioresult = zeros(numel(ratiofields), 3);
    
    for i = 1:numel(ratiofields)
        
       if positive(i) && design.(ratiofields{i}) < 0
           ratioresult(i,1) = true;
       end
           
       if zerotoone(i) && (design.(ratiofields{i}) < 0 || design.(ratiofields{i}) > 1)
           ratioresult(i,2) = true;
       end
       
       if ~isempty(ratiolimits)
          if (design.(ratiofields{i}) < ratiolimits(i,1) || design.(ratiofields{i}) > ratiolimits(i,2))
              ratioresult(i,3) = true;
          end
       end
        
    end

end

function s = keepfield(s, fields)
% removes all but the specified field(s) from a structure
%
% Syntax
% 
% s = keepfield(s, 'fieldname')
% s = keepfield(s, fields)
%
% Description
% 
% s = keepfield(s, 'fieldname') keeps the specified field in the structure
% array s, removing all others.
% 
% s = keepfield(s, fields) keeps several fields discarding all others. fields
% is a character array of field names or cell array of strings.
%
% See Also, FIELDNAMES, SETFIELD, GETFIELD, ISFIELD, ORDERFIELDS, RMFIELD

    if ischar(fields)
        fields = {fields};
    end
    
    
    if all(isfield(s, fields))
        
        fnames = fieldnames(s);
        
        for i = 1:numel(fields)
            
            fnames(strcmp(fnames, fields{i})) = [];
            
        end
        
        s = rmfield(s, fnames);
        
    else
        
        missingfields = fields(~isfield(s, fields));
        
        outstr = missingfields{1};
        for i = 2:numel(missingfields)
            outstr = [outstr, ', ', missingfields{i} ];
        end
        
        error('CROZIER:keepfields:InvalidFieldname', 'Field(s) %s not found in structure', outstr)
        
    end

end
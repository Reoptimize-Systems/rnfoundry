function ok = structHasAllFields (S, req_fields, throw, inputname)
% checks if all of a set of field names are present in a structure
%
% Syntax
%
%  ok = structHasAllFields (S, req_fields, throw)
%  ok = structHasAllFields (..., inputname)
%
% Input
%
%  S - structure to be tested for presence of fields. It will be checked if
%    all the fields in req_fields are present.
%
%  req_fields - cell string array of strings containing the names of fields
%    which must be present in the stucture.
%
%  throw - logical flag determining whether an error is thrown by
%    checkOrientationDescription if input fails check
%
%  inputname - (optional) string with name to use in error thrown by
%    structHasAllFields (if 'throw' is true). The error will be "<name> was
%    missing the following required fields: 'missingField1',
%    'missingField2', 'missingField3'".
%
% Output
%
%  ok - logical flag indicating if check was passed
%

    if nargin < 4
        inputname = 'input';
    end
    
    assert (isstruct (S), 'S must be a structure');
    if ~iscellstr (req_fields)
        error ('req_fields must be a cell array of strings containing a list of required field names in the input structure');
    end

    ok = true;
    
    tf = isfield (S, req_fields);
    
    if ~all(tf)
        ok = false;
    end

    if throw && ~ok
        
        missing_fields = req_fields(~tf);
        error ('%s was missing the following required fields: %s ', inputname, sprintf (['''%s''', repmat(', ''%s''', 1, numel(missing_fields)-1)], missing_fields{:}));

    end

end
function s = setfieldifabsent(s, field, v)
% sets the value of a structure field if, and only if, the field is not
% already present in the structure
%
% Syntax
%
% s = setfieldifabsent(s, field, v)
%
% Input
%
% s - scalar structure
%
% field - string containing the name of the structure field to be set to
%   the supplied value
%
% v - new contents of the structure field
%
% Output
%
% s - a structure with the field supplied in 'field' having contents 'v',
%   if it was not already present in the structure
%

    if ~ischar(field)
        error('field must be a character array');
    end
    
    if ~isfield(s, field)
        
        s.(field) = v;
        
    end

end
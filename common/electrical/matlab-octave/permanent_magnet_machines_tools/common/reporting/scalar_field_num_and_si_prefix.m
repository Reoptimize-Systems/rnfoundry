function [num, prefix] = scalar_field_num_and_si_prefix (S, fieldname, scalefac, idx)
% Converts a scalar numeric structure field into an SI prefixed string.
%
% Syntax
%
% [num, prefix] = scalar_field_num_and_si_prefix (S, fieldname)
% [...] = scalar_field_num_and_si_prefix (..., scalefac)
% [...] = scalar_field_num_and_si_prefix (..., scalefac, idx)
%
% Description
%
% scalar_field_num_and_si_prefix converts a scalar numeric structure field
% into an SI prefixed string. The value is shown in the string as a
% coefficient and an SI unit prefix, optimally chosen for readability. If
% the rounded |val|<10^-24 or |val|>=10^27 then E-notation is used, without
% a prefix. e.g. a force of 12000 Newtons stored in a field would be
% converted to '12 kN'.
%
% The following SI prefix string and symbols are used:
%
% Order  |1000^1 |1000^2 |1000^3 |1000^4 |1000^5 |1000^6 |1000^7 |1000^8 |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | kilo  | mega  | giga  | tera  | peta  | exa   | zetta | yotta |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|   k   |   M   |   G   |   T   |   P   |   E   |   Z   |   Y   |
%
% Order  |1000^-1|1000^-2|1000^-3|1000^-4|1000^-5|1000^-6|1000^-7|1000^-8|
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Name   | milli | micro | nano  | pico  | femto | atto  | zepto | yocto |
% -------|-------|-------|-------|-------|-------|-------|-------|-------|
% Symbol*|   m   |   u   |   n   |   p   |   f   |   a   |   z   |   y   |
%
% Input
%
%  S - structure containing field to be converted to a nicely formatted
%    number using SI prefixes for the power. (see sipre.m for more
%    information on the conversion).
%
%  fieldname - character vector containing the name of the field containing
%    the number to be converted.
%
%  scalefac - (optional) scale factor to apply to the number before doing
%    the conversion. If not supplied, default is 1 (no scaling).
%
%  idx - (optional) the structure field may be a vector or matrix. By
%    default scalar_field_num_and_si_prefix always converts only the number
%    at index 1 in the numeric field. The idx option can be used to convert
%    a number at a different index. Default is 1.
%
% Output
%
%  num - The number converted to a char vector (without SI prefix). If the
%    field was not present or not numeric, this will be the string 'N\A'.
%
%  prefix - The appropriate SI prefix for the number. If the field was not
%    present or not numeric this will be an empty string.
%
%
% See Also: sipre.m
%

    if nargin < 3
        scalefac = 1;
    else
        check.isNumericScalar (scalefac, true, 'scalefac');
    end
    
    if nargin < 4
        idx = 1;
    else
        check.isScalarInteger (idx, true, 'idx');
    end
    
    assert (isstruct (S) && isscalar (S), 'S must be a scalar structure');
    assert (ischar (fieldname), 'fieldname must be a character vector');
    
    
    if isfield (S, fieldname) ...
            && isnumeric (S.(fieldname))
        
        [numstr, prefix] = sipre (S.(fieldname)(idx) * scalefac, 3);
        
        num = numstr(1:end-numel(prefix)-1);
        
    else
        num = 'N/A';
        prefix = '';
    end

end
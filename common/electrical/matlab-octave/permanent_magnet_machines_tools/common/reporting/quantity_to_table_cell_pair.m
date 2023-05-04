function tcells = quantity_to_table_cell_pair (name, unit, S, fieldname, scalefac, idx)
% creates well formatted name-value table enrtries from structure fields
%
% Syntax
%
% tcells = quantity_to_table_cell_pair (name, unit, S, fieldname, scalefac)
%
% Description
%
% created a two element cell array with a name and value with an
% appropriate SI prefix, extracting the value from the field of a
% structure. Typically for use in conjuction with 
%
% Input
%
%  name - name of the quantity being entered in the table
%
%  unit - SI unit of the quantity
%
%  S - structure contining the quantity as a field.
%
%  fieldname - field name of the quantity in the structure, S
%
%  scalefac - (optional) scale factor to apply to the number before doing
%    the conversion. If not supplied, default is 1 (no scaling).
%
%  idx - (optional) the structure field may be a vector or matrix. By
%    default quantity_to_table_cell_pair always converts only the number
%    at index 1 in the numeric field. The idx option can be used to convert
%    a number at a different index. Default is 1.
%
%
% Output
%
%  tcells - a (1 x 2) cell array. The fisrt cell will contain the name of
%   the quantity (simply copied from the input), the second will contain an
%   SI prefixed string as produced by the function
%   "scalar_field_num_and_si_prefix".
%
%
%
% See Also: scalar_field_num_and_si_prefix.m
%

    if nargin < 5
        scalefac = 1;
    end

    if nargin < 6
        idx = 1;
    end

    [num, prefix] = scalar_field_num_and_si_prefix (S, fieldname, scalefac, idx);

    switch prefix

        case 'y'

            prefix = 'y';
        case 'z'
            prefix = 'z';
        case 'a'
            prefix = 'a';
        case 'f'
            prefix = 'f';
        case 'p'
            prefix = 'p';
        case 'n'
            prefix = 'n';
        case 'u'
            prefix = '$\mu$';
        case 'm'
            prefix = 'm';
        case 'k'
            prefix = 'k';
        case 'M'
            prefix = 'M';
        case 'G'
            prefix = 'G';
        case 'T'
            prefix = 'T';
        case 'P'
            prefix = 'P';
        case 'E'
            prefix = 'E';
        case 'Z'
            prefix = 'Z';
        case 'Y'
            prefix = 'Y';
    end

    if strcmp ('N/A', num)
        tcells = {name, num};
    else
        tcells = {[name, ' (', prefix, unit, ')'], num };
    end

end
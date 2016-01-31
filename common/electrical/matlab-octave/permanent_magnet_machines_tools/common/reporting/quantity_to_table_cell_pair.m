function tcells = quantity_to_table_cell_pair (name, unit, S, fieldname, scalefac)
% created a two element cell array with a name and value with an
% appropriate prefix


    if nargin < 5
        scalefac = 1;
    end

    [num, prefix] = scalar_field_num_and_si_prefix (S, fieldname, scalefac, 1);

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
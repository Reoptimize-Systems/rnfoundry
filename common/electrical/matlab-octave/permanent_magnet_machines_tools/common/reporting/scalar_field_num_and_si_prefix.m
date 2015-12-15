function [num, prefix] = scalar_field_num_and_si_prefix (S, fieldname, scalefac, idx)

    if nargin < 3
        scalefac = 1;
    end
    
    if nargin < 4
        idx = 1;
    end
    
    if isfield (S, fieldname) ...
            && isnumeric (S.(fieldname))
        
        [numstr, prefix] = sipre (S.(fieldname)(idx) * scalefac,3);
        
        num = numstr(1:end-numel(prefix)-1);
        
    else
        num = 'N/A';
        prefix = '';
    end

end
function str = floatforfilename (x)
    
    str = formatNumber (x);
    
    str = strrep (str, '.', 'pt');
    
end

function numstr = formatNumber (num)
% fomats a decimal number for pretty output to mbdyn file
%
% If it is an integer (to machine precision), it will be output
% with one decimal place. If a float it will be output to 14
% decimal places, but with any trailing zeros stripped to
% shorten it.

    if check.isint2eps (num)
        numstr = sprintf ('%d.0', num);
    else
        numstr = sprintf ('%.18f', num);
        % strip trailing zeros from decimals
        n = numel (numstr);
        while numstr(n) ~= '.'
            if numstr(n) == '0'
                numstr(n) = [];
            else
                break;
            end
            n = n - 1;
        end
    end

end
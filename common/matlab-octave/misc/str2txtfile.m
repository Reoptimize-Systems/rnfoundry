function str2txtfile (file, str)
% Write a string or character vector to a text file
%
% Syntax
%
% str2txtfile (file, strs)
%

    if isstring (str)
        cellstr2txtfile (file, str);
    elseif iscellstr (str)
        cellstr2txtfile (file, str);
    else
        cellstr2txtfile (file, {str});
    end
    
end
function range = xlsrange(row1, col1, row2, col2)
% constucts a range in excel syntax from row and column numbers
%
% Syntax
%
% range = xlsrange(row, col)
% range = xlsrange(row1, col1, row2, col2)
%
% Description
%
% xlsrange generates cell positions and ranges in excel's letter-number
% format, e.g. A1, B5:D9 etc. from row and column numbers.
%
% Examples
%
% range = xlsrange(1, 1)
% 
% range =
% 
% A1
%
% 
% range = xlsrange(1, 1, 4, 3)
% 
% range =
% 
% A1:C4
%

% Copyright Richard Crozier 2013


    if nargin == 2
        
        if ~(isscalar(row1) && isscalar(col1) && isnumeric(row1) && isnumeric(col1))
            error('xlsrange only supports numeric scalar inputs.')
        end
        
        col1letters = xlscol(col1);
        row1str = num2str(row1);
        
        range = [col1letters, row1str];
        
    elseif nargin == 4
        
        range = [ xlsrange(row1, col1), ':' xlsrange(row2, col2) ];
        
    else
        
        error('You must supply two or four arguments.');
        
    end
    
end
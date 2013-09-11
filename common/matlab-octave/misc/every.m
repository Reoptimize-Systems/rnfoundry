function tf = every(A)
% returns true if every element in a matrix evaluates to true
%
% See also: all

% Created by Richard Crozier 2011

    if isnumeric(A) || islogical(A)
        
        tf = all(A(:));
    
    else
        error('A must be a numeric or logical matrix');
    end
    
end
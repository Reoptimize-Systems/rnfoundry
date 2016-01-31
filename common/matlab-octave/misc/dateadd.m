function R = dateadd(D, V)
% dateadd: adds a period of time to a date
%
% Syntax
% 
% R = dateadd(D, V)
% 
% Input
%
% D - a date, in any format recognised by matlab
%
% V - a six element vector containing the amounts of time to be added where
%     
%     V  = [years, months, days, hours, minutes, seconds]
%
% Output
%
% R - A serial date number representing the supplied date pluss the
%     additional time
%
    
    if ~isscalar (D)
        R = datenum (D);
    else
        R = D;
    end
    
    R = addtodate (R, V(1), 'year');
    
    R = addtodate (R, V(2), 'month');
    
    R = addtodate (R, V(3), 'day');

    R = R + datenum (0,0,0,V(4),V(5),V(6)) - datenum (0,0,0,0,0,0);

end
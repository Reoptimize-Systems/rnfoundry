function result = evalhppolyeven(p, xp, halfperiod, xo)
% evaluates an even periodic function given the size of half a period and a
% polynomial fitted to the data from one half period of the function under
% the assumption that the function is inverted about the half period point,
% e.g. sin(x)
%

    if nargin < 4
        xo = [];
    end

    remxVhalfperiod = rem((xp./halfperiod),1);
    
    polysign = ones(size(xp,1),1);

    polysign(round( rem( (xp/halfperiod) - remxVhalfperiod, 2) ) == 0) = -1;
    
%     % check if period number is even or odd
%     if round( rem( (xp/halfperiod) - remxVhalfperiod, 2) ) == 0
%         polysign = -1;
%     else
%         polysign = 1;
%     end

    

    if polysign == 1
        
        if xp < 0
            result = polyvaln(p, [1 - abs(rem(xp,halfperiod)./halfperiod), xo]);
        else
            result = -polyvaln(p, [1 - abs(rem(xp,halfperiod)./halfperiod), xo]);
        end
        
    else
        
        if xp < 0
            result = -polyvaln(p, [abs(rem(xp,halfperiod)./halfperiod), xo]);
        else
            result = polyvaln(p, [abs(rem(xp,halfperiod)./halfperiod), xo]);
        end
        
    end
    
end
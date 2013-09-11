function result = evalhpslmodd(slm, evalmode, xp, halfperiod)
% evaluates an odd periodic function given the size of half a period and a
% polynomial fitted to the data from one half period of the function under
% the assumption that the function is inverted about the half period point,
% e.g. sin(x)
%

    remxVhalfperiod = rem((xp./halfperiod),1);
    
    polysign = ones(size(xp,1),1);

    polysign(round( rem( (xp/halfperiod) - remxVhalfperiod, 2) ) == 0) = -1;
    
    xR = abs(rem(xp,halfperiod)./halfperiod);
    
    xR(polysign == 1) = 1 - xR(polysign == 1);
    
    result = zeros(size(xp));
    
    if sum(xp<0) > 0
        result(xp<0) = polysign(xp<0,:) .* slmeval(xR(xp<0, :), slm, evalmode, false);

    end

    if sum(xp>=0) > 0
        result(xp>=0) = -polysign(xp>=0,:) .* slmeval(xR(xp>=0, :), slm, evalmode, false);
    end
    
    
%     % check if period number is even or odd
%     if round( rem( (xp/halfperiod) - remxVhalfperiod, 2) ) == 0
%         polysign = -1;
%     else
%         polysign = 1;
%     end
%
%     if polysign == 1
%         
%         if xp < 0
%             result = polyvaln(p, [1 - abs(rem(xp,halfperiod)./halfperiod), xo]);
%         else
%             result = -polyvaln(p, [1 - abs(rem(xp,halfperiod)./halfperiod), xo]);
%         end
%         
%     else
%         
%         if xp < 0
%             result = -polyvaln(p, [abs(rem(xp,halfperiod)./halfperiod), xo]);
%         else
%             result = polyvaln(p, [abs(rem(xp,halfperiod)./halfperiod), xo]);
%         end
%         
%     end
    
end
function tf = ismonatonicvec (V, increasing)
% checks if a vector is monatonically increasing or decreasing

    if nargin < 2
        increasing = true;
    end
    
    if isempty (V) || isscalar (V)
        tf = false;
        return;
    end
    
    if ~isvector (V)
       error ('Input is not a vector.') ;
    end
    
    temp = V(2:end) - V(1:end-1);
    
    if increasing
        if any (temp <= 0)
            tf = false;
            return;
        end
    else
        if any (temp >= 0)
            tf = false;
            return;
        end
    end

    tf = true;

end

function [pos, delt] = GetLinearMotionPos(Taup, v, n, minp, mint)
    
    tDist = minp * Taup;
    
    t = tDist / v;
    
    if t < mint
        tDist = v * mint;
    end
    
    dx = Taup / n;
    
    delt = dx / v;
    
    pos = 0:dx:tDist;

end
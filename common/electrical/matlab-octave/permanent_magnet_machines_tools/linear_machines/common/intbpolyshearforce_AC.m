function force = intbpolyshearforce_AC(design, J, xtemp)

    force = zeros(size(xtemp));

    x = xtemp(round(rem(xtemp-rem(xtemp, 1), 2)) == 0 & xtemp > 0);

    if ~isempty(x)
        force(round(rem(xtemp-rem(xtemp,1), 2)) == 0 & xtemp > 0) = ylorentzforce(J(round(rem(xtemp-rem(xtemp, 1), 2)) == 0 & xtemp > 0), design.p_intBx, abs(rem(x, 1)), design.MTL);
    end

    x = xtemp(round(rem(xtemp-rem(xtemp,1), 2)) ~= 0 & xtemp > 0);

    if ~isempty(x)
        force(round(rem(xtemp-rem(xtemp,1), 2)) ~= 0 & xtemp > 0) = -ylorentzforce(J(round(rem(xtemp-rem(xtemp,1), 2)) ~= 0 & xtemp > 0), design.p_intBx, abs(rem(x, 1)), design.MTL);
    end
    
    x = xtemp(round(rem(xtemp-rem(xtemp, 1), 2)) == 0 & xtemp < 0);

    if ~isempty(x)
        force(round(rem(xtemp-rem(xtemp,1), 2)) == 0 & xtemp < 0) = -ylorentzforce(J(round(rem(xtemp-rem(xtemp, 1), 2)) == 0 & xtemp < 0), design.p_intBx, abs(rem(x, 1)), design.MTL);
    end

    x = xtemp(round(rem(xtemp-rem(xtemp,1), 2)) ~= 0 & xtemp < 0);

    if ~isempty(x)
        force(round(rem(xtemp-rem(xtemp,1), 2)) ~= 0 & xtemp < 0) = ylorentzforce(J(round(rem(xtemp-rem(xtemp,1), 2)) ~= 0 & xtemp < 0), design.p_intBx, abs(rem(x, 1)), design.MTL);
    end

    
end
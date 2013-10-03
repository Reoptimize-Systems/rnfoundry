function psidot = polypsidot_ACTM(design, xR)

    xtemp = xR;

    psidot = zeros(size(xtemp));

    x = xtemp(round(rem((xtemp)-rem((xtemp), 1), 2)) == 0 & xtemp > 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) == 0 & xtemp > 0) = polyvaln(design.p_psidot, abs(rem(x, 1)));
    end

    x = xtemp(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp > 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp > 0) = -polyvaln(design.p_psidot, abs(rem(x, 1)));
    end
    
    x = xtemp(round(rem((xtemp)-rem((xtemp), 1), 2)) == 0 & xtemp < 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) == 0 & xtemp < 0) = -polyvaln(design.p_psidot, abs(rem(x, 1)));
    end

    x = xtemp(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp < 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp < 0) = polyvaln(design.p_psidot, abs(rem(x, 1)));
    end
    
end
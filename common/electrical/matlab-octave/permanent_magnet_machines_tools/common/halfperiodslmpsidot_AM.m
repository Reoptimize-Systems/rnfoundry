function psidot = halfperiodslmpsidot_AM(design, xR, xscale)
% calculates the derivative of the coil flux linkage defined by an slm
% object fitted to a half period of the flux linkage waveform
%

    xtemp = xR;

    psidot = zeros(size(xtemp));

    x = xtemp(round(rem((xtemp)-rem((xtemp), 1), 2)) == 0 & xtemp > 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) == 0 & xtemp > 0) = slmeval(abs(rem(x, 1)), design.slm_psidot, 1, false);
    end

    x = xtemp(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp > 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp > 0) = -slmeval(abs(rem(x, 1)), design.slm_psidot, 1, false);
    end
    
    x = xtemp(round(rem((xtemp)-rem((xtemp), 1), 2)) == 0 & xtemp < 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) == 0 & xtemp < 0) = -slmeval(abs(rem(x, 1)), design.slm_psidot, 1, false);
    end

    x = xtemp(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp < 0);

    if ~isempty(x)
        psidot(round(rem((xtemp)-rem((xtemp),1), 2)) ~= 0 & xtemp < 0) = slmeval(abs(rem(x, 1)), design.slm_psidot, 1, false);
    end
    
    psidot = psidot ./ xscale;
    
end
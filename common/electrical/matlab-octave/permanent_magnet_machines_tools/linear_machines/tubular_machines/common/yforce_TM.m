function Fy = yforce_TM(design, Jz, y)

    
    ytemp = y ./ design.Wp;
    
    Fy = zeros(size(ytemp));

    y = ytemp(round(rem((ytemp)-rem((ytemp), 1), 2)) == 0 & ytemp > 0);

    if ~isempty(y)
        Jztemp = Jz(round(rem((ytemp)-rem((ytemp), 1), 2)) == 0 & ytemp > 0);
        Fy(round(rem((ytemp)-rem((ytemp),1), 2)) == 0 & ytemp > 0) = ylorentzforce(Jztemp, design.p_intBx, abs(rem(y, 1)), design.MTL);
    end

    y = ytemp(round(rem((ytemp)-rem((ytemp),1), 2)) ~= 0 & ytemp > 0);

    if ~isempty(y)
        Jztemp = Jz(round(rem((ytemp)-rem((ytemp),1), 2)) ~= 0 & ytemp > 0);
        Fy(round(rem((ytemp)-rem((ytemp),1), 2)) ~= 0 & ytemp > 0) = -ylorentzforce(Jztemp, design.p_intBx, abs(rem(y, 1)), design.MTL);
    end
    
    y = ytemp(round(rem((ytemp)-rem((ytemp), 1), 2)) == 0 & ytemp < 0);

    if ~isempty(y)
        Jztemp = Jz(round(rem((ytemp)-rem((ytemp), 1), 2)) == 0 & ytemp < 0);
        Fy(round(rem((ytemp)-rem((ytemp),1), 2)) == 0 & ytemp < 0) = -ylorentzforce(Jztemp, design.p_intBx, abs(rem(y, 1)), design.MTL);
    end

    y = ytemp(round(rem((ytemp)-rem((ytemp),1), 2)) ~= 0 & ytemp < 0);

    if ~isempty(y)
        Jztemp = Jz(round(rem((ytemp)-rem((ytemp),1), 2)) ~= 0 & ytemp < 0);
        Fy(round(rem((ytemp)-rem((ytemp),1), 2)) ~= 0 & ytemp < 0) = ylorentzforce(Jztemp, design.p_intBx, abs(rem(y, 1)), design.MTL);
    end
    
    %Fy = sum(ylorentzforce(Jz, design.p_intBx, y));
    
    %Fy = Fy .* design.Phases;
    
    Fy = sum(Fy,1) .* design.Poles(1);

end
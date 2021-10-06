function psi = psilookuppsi_linear(design, y)

    psi = zeros(size(y));

    for i = 1:length(y)

        if (y(i) < 0 && round(rem(y(i)-rem(y(i), 1), 2)) == 0) || (y(i) > 0 && ~(round(rem(y(i)-rem(y(i), 1), 2)) == 0))
            psi(i) = -interp1(design.psilookup(1,:), design.psilookup(2,:), abs(rem(y(i), 1)));
        else
            psi(i) = interp1(design.psilookup(1,:), design.psilookup(2,:), abs(rem(y(i), 1)));
        end

    end

end
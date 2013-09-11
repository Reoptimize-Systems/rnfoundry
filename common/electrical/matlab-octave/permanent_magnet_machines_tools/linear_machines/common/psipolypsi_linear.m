function psi = psipolypsi_linear(design, y)

    psi = zeros(size(y));

    for i = 1:length(y)

        if round(rem((y(i))-rem((y(i)), 1), 2)) == 0
            psi(i) = polyvaln(design.psipoly, abs(rem(y(i), 1)));
        else
            psi(i) = polyvaln(design.psipoly, 1-abs(rem(y(i), 1)));
        end

    end

end
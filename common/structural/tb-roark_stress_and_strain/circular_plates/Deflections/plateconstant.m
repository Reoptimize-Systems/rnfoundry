function D = plateconstant(E, t, v)

    D = E .* t.^2 ./ (12 .* (1 - v.^2));

end
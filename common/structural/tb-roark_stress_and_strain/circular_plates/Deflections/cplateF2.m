function F2 = cplateF2(b,r)

    F2 = 0.25 .* (1 - (b./r).^2 .* (1 + 2.*log(r./b)));

end
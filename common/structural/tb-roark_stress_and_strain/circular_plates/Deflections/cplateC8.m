function c8 = cplateC8(a, b, v)

    c8 = 0.5 .* ( 1 + v + (1 - v).*(b ./ a).^2 );

end
function c2 = cplateC2(a, b)

    c2 = 0.25 .* ( 1 - (b./a).^2 .*(1 + 2.*log(a./b)) );

end
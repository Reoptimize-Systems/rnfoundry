function c9 = cplateC9(a, b, v)

    c9 = (b./a) .* ( ((1+v)./2) .* log(a./b) + ((1-v)./4).*(1 - (b./a).^2) );

end
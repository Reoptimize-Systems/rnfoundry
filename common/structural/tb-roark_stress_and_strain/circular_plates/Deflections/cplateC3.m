function c3 = cplateC3(a, b)

    c3 = (b ./ (4.*a)) .* ( ((b./a).^2 + 1).*log(a./b) + (b./a).^2 - 1);

end
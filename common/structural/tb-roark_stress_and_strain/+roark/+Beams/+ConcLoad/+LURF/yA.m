function res = yA (W,l,a,E,I)

    res = -W .* ( 2.*realpow (l,3) - 3.*realpow (l,2).*a + realpow (a,3) ) ./ (6.*E.*I);
    
end
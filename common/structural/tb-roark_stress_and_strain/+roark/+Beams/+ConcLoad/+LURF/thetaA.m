function res = thetaA (W,l,a,E,I)

    res = W .* realpow(l-a, 2) ./ (2.*E.*I) ;

end
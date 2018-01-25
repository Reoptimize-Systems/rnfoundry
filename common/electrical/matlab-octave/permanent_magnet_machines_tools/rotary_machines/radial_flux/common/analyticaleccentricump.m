function Fnet = analyticaleccentricump (g, R, L, dg, Bipeak)
% analytically calculates the unbalanced magnet pull for an eccentric rotor
    
    epsilon = dg ./ g;

    Fnet = realpow (Bipeak, 2) .* epsilon .* R .* L .* pi ./ (2 .* mu_0);
    
end
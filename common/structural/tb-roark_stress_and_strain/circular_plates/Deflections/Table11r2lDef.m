function y = Table11r2lDef(Gvars, Lvars, E, v, r)

    a = Gvars(:,1);
    b = Gvars(:,2);
    t = Gvars(:,3);
    
    r_o = Lvars(:,1);
    q = Lvars(:,2);
    
    D = plateconstant(E, t, v);
    
    yb = 0;
    
    thetab = 0;
    
    Mrb = (-q .* a.^2 ./ cplateC8(a, b, v)) ...
           .* (((cplateC9(a, b, v) ./ (2.*a.*b)).*(a.^2 - r_o.^2)) - cplateL17(r_o, a, v));
    
    Qb = (q ./ (2.*b)) .* (a.^2 - r_o.^2);

%     ya = Mrb .* cplateC2(a, b) .* (realpow(a,2)./D) ...
%         + Qb .* cplateC3(a, b) .* (a.^3 ./ D) ...
%         - cplateL11(r_o, a) .* (q .* a.^4 ./ D);
    
    F1 = cplateF1(a,b,v,r);
    
    F2 = cplateF2(b,r);
    
    F3 = cplateF3(b,r);

    y = zeros(size(Gvars,1), numel(r));

    for i = 1:size(Mrb, 1)
        y(i,:) = GenericYDefAnnularPlateUniformDistPressure( ...
                yb, D, thetab, F1(i), Mrb(i), F2(i), Qb(i), F3(i), q(i), cplateG11(r_o(i),r), r);
    end

end





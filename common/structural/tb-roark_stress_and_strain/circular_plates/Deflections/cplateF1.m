function F1 = cplateF1(a,b,v,r)

    F1 = ((1 + v.*a.*b) ./ 2).*(b./r).*log(r./b) + ((1-v)./4).*((r./b)-(b./r));

end
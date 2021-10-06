function F3 = cplateF3(b,r)

    bVr = b./r;
    
    F3 = (b./(4.*r)) .* ((bVr.^2 + 1).*log(r./b) + (bVr).^2 - 1);
    
end


function Ewtot = GetAnnualEnergy(Betah_data, Betag, TpVHs_fdata, TpVHs_fscale, Hs, Tw)

    TpVHs_fdata = TpVHs_fdata .* TpVHs_scale;
    
    Hs = repmat(Hs', 1, size(Betah_fdata, 2));
    
    Tw = repmat(Tw, size(Betah_fdata, 1), 1);
    
    Ewtot = ((Betah_data.^2) + (Betah_data.*Betag) + (Betag.^2)) .* 2 .* (pi^2) .* (Hs.^2) ./ ((Betah_data + Betag + 1) .* Tw);
    
    Ewtot = Ewtot .* TpVHs_fdata .* TpVHs_fscale;
    
    Ewtot = sum(Ew);
    
    Ewtot = sum(Ew);
    
end
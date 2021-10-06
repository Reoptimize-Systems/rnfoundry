function psi = fluxlinkagefromAdata_TM(design, x)
% calculates the flux linkage in a coil of a tubular machine from sampled
% vector potential data
%
% Syntax
%
% psi = fluxlinkagefromAdata_TM(design, x)
%
% 

    ymin = x - (design.WcVWp / 2);
    ymax = x + (design.WcVWp / 2);

    intA = zeros(size(x));
    
    [intA(1), intAslm] = integratehalfperiod2ddata( design.X, design.Y, design.A, ...
                                                    design.Ri, ...
                                                    ymin(1), ...
                                                    design.Ro - design.FEMMTol - 2*max(eps(design.X(:))), ...
                                                    ymax(1) );
    
    for i = 1:length(x)
        
         intA(i) = integratehalfperiod2ddata(intAslm, ymin(i), ymax(i));
        
    end
    
    psi = design.CoilTurns .* intA .* design.Wp ./ design.CoilArea; 
    
end
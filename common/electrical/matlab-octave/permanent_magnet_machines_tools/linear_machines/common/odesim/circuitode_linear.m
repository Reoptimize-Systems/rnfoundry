function dIphases = circuitode_linear(Iphases, EMF, design)
% circuitode_linear: solves the rhs of the differential equation describing
% the linear machine and grid circuit

    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
%     dIphases = ((design.CoilsPerBranch * EMF) - (Iphases .* design.R)) ./ design.L;
    
    dIphases = multiphasecurrentderiv(Iphases, EMF, design.RPhase + design.RLoad, design.L);
    
end
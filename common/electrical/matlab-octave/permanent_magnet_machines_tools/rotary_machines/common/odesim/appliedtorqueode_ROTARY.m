function [dx, varargout] = appliedtorqueode_ROTARY(t, x, design, simoptions, torque)

    if nargin == 0
        % get the output names
        varargout{1} = {'ydot', 'dpsidxR', 'EMF', 'xT', 'vT', 'Fpto', 'FaddE', 'FaddEBD', 'RPhase'};
        return;
    end
    
    xR = x(1);
    vR = x(2);
    
    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(3:(2+design.Phases));
    
    Icoils = Iphases ./ design.Branches;
    
    % preallocate dx
    dx = zeros(size(x));
    
    % determine the machine outputs
    [dpsidxR, EMF, Feff, FfeaVec, xT, vT, unitv, design] = ...
        machineodesim_linear(design, simoptions, Icoils, xR, vR, 0, 0);
    
    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(3:(2+design.Phases),1) = circuitode_linear(Iphases, EMF, design);
     
    % call the supplied additional force function
    [FaddE, ForceBD] = ...
        feval(simoptions.ODESim.TorqueFcn, design, simoptions, xT, vT, EMF, Iphases, simoptions.ODESim.TorqueFcnArgs{:});
    
    % velocity
    dx(1) = vR;
    
    % convert torque to tangential force
    Finput = torque / design.Rmm;
    
    % acceleration
    dx(2) = (Finput + Feff + FaddE) / design.RotorMass;
    
    % To record the forces
    if nargout > 1
        
        % rate of change of flux linkage with displacement 
        varargout{2} = dpsidxR;
        % per-coil induced voltage 
        varargout{3} = EMF;
        % vertical displacement of the translator  
        varargout{4} = xT;
        % vertical velocity of the translator
        varargout{5} = vT;
        % total electromagnetic forces acting on effector
        varargout{6} = Feff;
        % output the additional forces transmitted to the primemover
        varargout{7} = FaddE;
        % output a breakdown of the additional forces transmitted to the
        % prime mover
        varargout{8} = ForceBD;
        % output the phase resistance at each time step
        varargout{9} = diag(design.RPhase)';
        
    end

end
function varargout = feaprescribedmotodetorquefcn_ROTARY(t, x, design, simoptions)
% solves the right had side of the differential equations for a machine
% moving with a prescribed motion using FEA to determine main quantities
%
% Syntax:
%
% Input
%
% x is a vector of values of the state of the sytem at time t. The members
% of x are the following:
%
%
% design is a standard machine design structure popluated with all the
% necessary information to perform the simulation. 
%
% 'simoptions'  is a structure determining the buoy and sea parameters used
% in the simulation, and also any penalties to be used in the scoring of
% the design. Essential fields are:

    if nargin == 0
        % get the output names
        varargout{1} = {'ydot', 'dlambdaVdt', 'EMF', 'thetaT', 'omegaT', 'Tqpto', 'TqaddE', 'TqaddEBD', 'RPhase'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices,:);
    
    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));
    
    % get the velocity and position at the current time
    [thetaE, omegaE] = prescribedmotomegatheta(t, simoptions);

    % determine the machine outputs Get the flux linkage and forces using
    % the core machine simulation function
    [flux_linkage, Tqeff, Tqreac, lossinfo, design] = ...
        machinefeaodesim_AM (design, simoptions, thetaE, 0, omegaE, 0, Icoils);
    
    % clculate the numerical derivative of the flux linkage
    flodederiv = simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.NumericalDerivatives;
    dlambdaVdt = flodederiv.derivative (t, flux_linkage);
    
    % the emf is the negative of the flux linkage derivative
    EMF = -dlambdaVdt .* design.CoilsPerBranch;
    
    % do the integration
    dx(simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.SolutionIndices,:) = ...
        dlambdaVdt;
    
    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    for ind = 1:size (Iphases, 2)
        dx(simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices,ind) = ...
            circuitode_linear (Iphases(:,ind), EMF(:,ind), design);
    end
     
    % call the supplied additional force function
    for ind = 1:size (Iphases, 2)
        [TqaddE(ind), TorqueBD(ind,:)] = ...
            feval (simoptions.ODESim.TorqueFcn, design, simoptions, thetaE, omegaE, EMF(:,ind), Iphases(:,ind), simoptions.ODESim.TorqueFcnArgs{:});
    end
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % rate of change of flux linkage with displacement 
        varargout{2} = dlambdaVdt;
        % per-coil induced voltage 
        varargout{3} = EMF;
        % vertical displacement of the translator  
        varargout{4} = thetaE;
        % vertical velocity of the translator
        varargout{5} = omegaE;
        % total electromagnetic forces acting on effector
        varargout{6} = Tqeff;
        % output the additional forces transmitted to the primemover
        varargout{7} = TqaddE;
        % output a breakdown of the additional forces transmitted to the
        % prime mover
        varargout{8} = TorqueBD;
        % output the phase resistance at each time step
        varargout{9} = diag(design.RPhase)';
        
        
    end

end


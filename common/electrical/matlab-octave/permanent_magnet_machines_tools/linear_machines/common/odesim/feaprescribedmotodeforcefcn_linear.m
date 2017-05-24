function varargout = feaprescribedmotodeforcefcn_linear (t, x, design, simoptions)
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
        varargout{1} = {'ydot', 'dlambdaVdt', 'EMF', 'xT', 'vT', 'Fpto', 'FaddE', 'FaddEBD', 'RPhase'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices,:);
% Iphases = zeros (size (x));
    
    Icoils = Iphases ./ design.Branches;

    ncols = size (Iphases, 2);
    
    % Initialize dx with zeros
    dx = zeros (size (x));
    
    % get the velocity and position at the current time
    [xE, vE] = prescribedmotvelpos (t, simoptions);

    % determine the machine outputs Get the flux linkage and forces using
    % the core machine simulation function
    [flux_linkage, Feff, Freac, lossinfo, design] = ...
        machinefeaodesim_AM (design, simoptions, xE, 0, vE, 0, Icoils);
    
    fprintf (1, 'tguess: %g, flux_linkage: %s, Icoils: %s\n', t, mat2str (flux_linkage'), mat2str (Icoils'));
    
%     [flux_linkage, Feff, Freac, lossinfo, design] = ...
%         slmflmachineodesim_AM (design, simoptions, ...
%                                repmat (xE, 1, ncols), 0, ...
%                                repmat (vE, 1, ncols), 0, ...
%                                Icoils );
    
%     if t == simoptions.ODESim.TimeSpan(1)
%         [~, ~, ~, dlambdaVdt, design, ~] = ...
%             machineodesim_AM ( design, simoptions, ...
%                                repmat (xE + design.zp/2 + 2*design.zp/3, 1, ncols), 0, ...
%                                repmat (vE, 1, ncols), 0, ...
%                                Icoils );
% %         dlambdaVdt = slmeval
%     else
        % clculate the numerical derivative of the flux linkage
        flodederiv = simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.NumericalDerivatives;
        dlambdaVdt = flodederiv.derivative (t, flux_linkage);
%     end
    
    % the emf is the negative of the flux linkage derivative
    EMF = -dlambdaVdt .* design.CoilsPerBranch; % * -1;
    
    % do the integration
    dx(simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.SolutionIndices,:) = ...
        dlambdaVdt;
    
%     % find the derivative of the coil current (solving the differential
%     % equation describing the simple output circuit)
    for ind = 1:size (Iphases, 2)
%         if all (EMF) == 0
%             dx(simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices,ind) = 0;
%         else
            dx(simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices,ind) = ...
                circuitode_linear (Iphases(:,ind), EMF(:,ind), design);
%         end
    end
     
    % call the supplied additional force function
    for ind = 1:size (Iphases, 2)
        [FaddE(ind), ForceBD(ind,:)] = ...
            feval (simoptions.ODESim.ForceFcn, design, simoptions, xE, vE, EMF(:,ind), Iphases(:,ind), simoptions.ODESim.ForceFcnArgs{:});
    end
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % rate of change of flux linkage with displacement 
        varargout{2} = dlambdaVdt;
        % per-phase induced voltage 
        varargout{3} = EMF;
        % vertical displacement of the translator  
        varargout{4} = xE;
        % vertical velocity of the translator
        varargout{5} = vE;
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


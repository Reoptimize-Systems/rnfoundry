function varargout = prescribedmotodeforcefcn_linear(t, x, design, simoptions)
% prescribedmotodeforcefcn_linear: solves the right had side of the differential
% equations for a machine moving with a prescribed motion
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
%
% TODO: complete help
%

    if nargin == 0
        % get the output names
        varargout{1} = {'ydot', 'dpsidxR', 'EMF', 'xT', 'vT', 'Fpto', 'FaddE', 'FaddEBD', 'RPhase'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(1:design.Phases);

    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));
    
    % get the velocity and position at the current time
    [xTtemp, vTtemp] = prescribedmotvelpos(t, simoptions);

    simoptions.tether_length = 1000;

    % determine the machine outputs
    [dpsidxR, EMF, Feff, FfeaVec, xT, vT, unitv, design] = ...
        machineodesim_linear(design, simoptions, Icoils, xTtemp, vTtemp, 0, 0);
    
    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(1:design.Phases,1) = circuitode_linear(Iphases, EMF, design);
     
    % call the supplied additional force function
    [FaddE, ForceBD] = ...
        feval(simoptions.forcefcn, design, simoptions, xT, vT, EMF, Iphases, simoptions.forcefcnargs{:});
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

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


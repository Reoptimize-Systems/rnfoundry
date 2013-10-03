function varargout = prescribedmotode_linear(t, x, design, simoptions)
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

    if nargin == 0
        % get the output names
        varargout{1} = {'ydot', 'dpsidxR', 'EMF', 'xT', 'vT', 'Fpto', 'RPhase'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(1:design.phases);
    
    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));
    
    % get the velocity and position at the current time
    [xTtemp, vTtemp] = prescribedmotvelpos(t, simoptions);
    
    simoptions.tether_length = 1000;

    % determine the machine outputs
    [dpsidxR, EMF, Fpto, FfeaVec, xT, vT, unitv, design] = ...
        machineodesim_linear(design, simoptions, Icoils, xTtemp, vTtemp, 0, 0);
    
    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(1:design.phases,1) = circuitode_linear(Iphases, EMF, design);
    
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
        varargout{6} = Fpto;
        % the phase resistance
        varargout{7} = diag(design.RPhase);
        
    end

end


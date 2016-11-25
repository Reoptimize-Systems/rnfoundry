function varargout = nestedodeforcefcn_linear (t, x, mc, flag, design, simoptions)
% nestedodeforcefcn_linear: solves the right had side of the differential
% equations for a machine moving with a prescribed motion
%
% Syntax:
%
% Input
%
%  x is a vector of values of the state of the sytem at time t
%
%  mc - coefficients of linear interpolation functions of the form 
%    y = m*t + c for the position and velocity at the current time provided
%    as a two column matrix where mc(1,1:2) are the m and c values for the
%    position interpolation function and mc(2,1:2) are the m and c values
%    for the velocity interpolation function.
%
%  flag - scalar integer flag determining solution stage and output from
%    the function. If flag == 0 a cell array of strings containing the
%    names of internally calculated values will be returned. These
%    variables are returned by the function if more than one output is
%    requested. The order of the cell array is the same as orer of retunred
%    arguments. This is intended to assist with regenerating the internally
%    calulated values after the main solution is complete.
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

    if flag == 0
        % get the output names
        varargout{1} = {'ydot', 'dpsidxR', 'EMF', 'Fpto', 'FaddE', 'FaddEBD', 'RPhase', 'xE', 'vE'};
        return;
    end

    % Change the x members into more useful variables names, MATLAB will
    % optimise away any memory penalty associated with this I think    
    Iphases = x(1:design.Phases);

    Icoils = Iphases ./ design.Branches;

    % Initialize dx with zeros
    dx = zeros(size(x));
    
    % get the velocity and position at the current time
    xT = mc(1,1) .* t + mc(1,2); 
    vT = mc(2,1) .* t + mc(2,2);

    simoptions.BuoySim.tether_length = 1000;

    % determine the machine outputs
    [Feff, ~, EMF, dpsidxR, design] = machineodesim_AM (design, simoptions, xT, 0, vT, 0, Icoils);
%     [dpsidxR, EMF, Feff, FfeaVec, xT, vT, unitv, design] = ...
%         machineodesim_linear (design, simoptions, Icoils, xTtemp, vTtemp, 0, 0);
    
    % find the derivative of the coil current (solving the differential
    % equation describing the simple output circuit)
    dx(1:design.Phases,1) = ((EMF(:) - (design.RPhase + design.RLoad) * Iphases(:))' / design.L)'; % circuitode_linear (Iphases, EMF, design);
     
    % call the supplied additional force function
    [FaddE, ForceBD] = ...
        feval ( simoptions.ODESim.ForceFcn, ...
                design, simoptions, ...
                xT, vT, EMF, Iphases, ...
                simoptions.ODESim.ForceFcnArgs{:} );
    
    % ************************************************************************

    % Now assign the outputs
    varargout{1} = dx;

    % To record the forces
    if nargout > 1
        
        % rate of change of flux linkage with displacement 
        varargout{2} = dpsidxR;
        % per-coil induced voltage 
        varargout{3} = EMF;
        % total electromagnetic forces acting on effector
        varargout{4} = Feff;
        % output the additional forces transmitted to the primemover
        varargout{5} = FaddE;
        % output a breakdown of the additional forces transmitted to the
        % prime mover
        varargout{6} = ForceBD;
        % output the phase resistance at each time step
        varargout{7} = diag(design.RPhase)';
        % translator position
        varargout{8} = xT;
        % translator velocity
        varargout{9} = vT;
        
    end

end


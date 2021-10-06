function varargout = simplelinearmachineode_proscribedmotion(t, y, design, simoptions)
% General purpose function for solving the differential equations
% associated with the linear machines where a series of time and positions 
% for the translator are supplied
%
% 
%   y = [ I ]

    if nargin == 0

        % get the output names
        varargout{1} = {'ydot', 'dpsidxF', 'EMF', 'xT', 'vT', 'Force'};
        return;
    end
    
    % the translator is driven by a very stiff drive, and the displacement
    % and velocity calculated by interpolation of the values provided
    % in simoptions
    
    % Interpolate the data set (times,xT) at current time
    xTtemp = interp1(simoptions.drivetimes, simoptions.xT, t); 
    
    % Interpolate the data set (times,vT) at current time
    vTtemp = interp1(simoptions.drivetimes, simoptions.vT, t); 

    simoptions.BuoySim.tether_length = 1000;
    
    [dpsidxF, EMF, Force] = ...
        machinesim_linear(design, simoptions, y, xTtemp, 0, vTtemp, 0);
    
    % find the derivative of the coil current
    Idot = (EMF - y .* design.R) ./ design.L;

    ydot = Idot;

    varargout{1} = ydot;

    if nargout > 1
        varargout{2} = dpsidxF;
        varargout{3} = EMF;
        varargout{4} = xTtemp;
        varargout{5} = vTtemp;
        varargout{6} = Force;
    end

end






function varargout = slink_TM(y, v, I, design, dpsidxfun, Fyfun)

    % We will get the numerical derivatives of the flux linkage in each of
    % the coils at the current position
    y = [y, y + (design.Wp / 3), y + (2 * design.Wp / 3)];
    
    dpsidx = [feval(dpsidxfun, design, y(1)), ...
        feval(dpsidxfun, design, y(2)),...
        feval(dpsidxfun, design, y(3))];
    
    % determine the emf in the coils
    EMF = -v .* dpsidx;

    % find the current density in the coil block by finding the total
    % current passing through a section and dividing by the cross-sectional
    % area
    Jz = I .* design.CoilTurns ./ design.CoilArea;
    
    % calculate the force from the machine
    Fy = feval(Fyfun, design, Jz, y);

    varargout{1} = EMF;
    varargout{2} = Fy;

end






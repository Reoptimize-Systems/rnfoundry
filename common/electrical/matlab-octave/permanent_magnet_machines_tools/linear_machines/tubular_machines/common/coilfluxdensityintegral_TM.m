function [intBx, intBy] = coilfluxdensityintegral_TM(design, y, mode)
% Finds the integral of the flux density over the coil cross-section in the
% axial and radial directions for a tubular machine at multiple positions
%
% Syntax
% 
% [intBx, intBy] = coilfluxdensityintegral_TM(design, y, mode)
% 
% Input
%
% design is a structure containing at least the following fields:
%
%   WcVWp - length of region of integration in y direction relative to the
%     total length of region of smapled data in y direction (expected to be
%     one half period of periodic function (periodic in y). I.e the ratio
%     of coil width to pole width in the tubular machine.
%
%   Ri - lower bound of region of integration in x direction. Inner coil
%     radius for tubular machine
%
%   Ro - upper bound of region of integration in x direction. Outer coil
%     radius for tubular machine
%
% design must also have either the fields 'X', 'Y', 'Bx', 'By', or 'p_Bx' 
% and 'p_By'. If the first set are present these are sampled B data in the
% x and y direction and the respective x and y coordinates of each data
% point. If the second set these are polynomials fitted to the same data.
%


    if nargin < 3
    
        if every(isfield(design, {'X', 'Y', 'Bx', 'By'}))
            
            mode = 'data';
            
        elseif every(isfield(design, {'p_Bx', 'p_By'}))
            
            mode = 'poly';
            
        else
            error('You must provide either two polynomials fitted to the B data or the data as field of the design structure')
        end
        
    end
    
    ymin = y - (design.WcVWp / 2);
    ymax = y + (design.WcVWp / 2);
    
    xmin = repmat(design.Ri, size(y));
    xmax = repmat(design.Ro, size(y))-design.FEMMTol;
    
    mintol = 1e-9;
    
    tol = [mintol, mintol];
    
    edgemode = [0, 1];
    
    intBx = zeros(1, length(y));
    intBy = zeros(1, length(y));
    
    switch mode
        
        case 'data'
                        
            for i = 1:length(y)

                if i == 1
                    
                    [intBx(i), Bxintslm] = integratehalfperiod2ddata(design.X, design.Y, design.Bx, xmin(i), ymin(i), xmax(i), ymax(i), edgemode(1));

                    [intBy(i), Byintslm] = integratehalfperiod2ddata(design.X, design.Y, design.By, xmin(i), ymin(i), xmax(i), ymax(i), edgemode(2));

                else
                    intBx(i) = integratehalfperiod2ddata(Bxintslm, ymin(i), ymax(i));

                    intBy(i) = integratehalfperiod2ddata(Byintslm, ymin(i), ymax(i));
                    
                end
                
                tol = abs([intBx(i) / 1000, intBy(i) / 1000]);

                tol(tol<mintol) = mintol;

            end

        case 'poly'

            for i = 1:length(y)

                intBx(i) = integratehalfperiodypoly(design.p_Bx, xmin(i), ymin(i), xmax(i), ymax(i), tol(1), edgemode(1));

                intBy(i) = integratehalfperiodypoly(design.p_By, xmin(i), ymin(i), xmax(i), ymax(i), tol(2), edgemode(2));

                tol = abs([intBx(i) / 1000, intBy(i) / 1000]);

                tol(tol<mintol) = mintol;

            end

    end

    % data and polynomial are given for a y range between 0 and one, we
    % must scale the resulting integral by the pole width to denormalise
    intBx = intBx .* design.PoleWidth;

    intBy = intBy .* design.PoleWidth;
    
end
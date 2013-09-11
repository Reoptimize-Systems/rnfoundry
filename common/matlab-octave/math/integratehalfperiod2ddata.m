function [int, ylineintslm] = integratehalfperiod2ddata(X, Y, Z, xmin, ymin, xmax, ymax, mode)
% finds the double integral of a set of 2D data
%
% Syntax
%
% [int, ylineintslm] = integratehalfperiod2ddata(X, Y, Z, xmin, ymin, xmax, ymax, tol, mode)
% [int, ylineintslm] = integratehalfperiod2ddata(X, Y, Z, xmin, ymin, xmax, ymax, tol, mode, ylineintslm)
% [int, ylineintslm] = integratehalfperiod2ddata(ylineintslm, ymin, ymax)
%
% Input
%
%   X, Y, Z - nterpolation within the two-dimensional function specified by
%             matrices X, Y, and Z. X and Y must be monotonic, and have the
%             same format ("plaid") as if they were produced by meshgrid.
%             Matrices X and Y specify the points at which the data Z is
%             given
%
%   The values xmin, ymin, xmax, ymax describe a rectangular region
%   within which the integral of the polynomial is to be evaluated
%
%   xmin - minimum value of x
%
%   ymin - normalized maximum value of y so that the range 0 <= y < 1 represents
%          one half period of the actual function 
%
%   xmax - maximum value of x
%
%   ymax - normalized minimum value of y so that the range 0 <= y < 1 represents
%          one half period of the actual function 
%
%   tol
%
%   mode - (optional) provides information as to the type of function to be
%          integrated, 
    
    if nargin < 8
        mode = 1;
    end
    
    if nargin ~= 3
        
        if ischar(mode)
            
            if strcmp(mode, 'max')
                mode = true;
            elseif strcmp(mode, 'zero')
                mode = false;
            else
                error('Unrecognised mode string.')
            end
            
        end
        
        % first generate line integrals in the non-periodic direction for a
        % full period in the periodic direction
        nxinterppoints = 2*size(X,1);
        xi = linspace(xmin, xmax, nxinterppoints);
        nlineintpnts = 100;
        y = linspace(min(Y(:)), max(Y(:)) + (max(Y(:)) - min(Y(:))), nlineintpnts);
        for i = 1:nlineintpnts
%             if mode
%                 lineints(i) = quad(@(x) maxedgesinterpsubfcn(x, y(i), X, Y, Z), xmin, xmax, tol);
%             else
%                 lineints(i) = quad(@(x) zeroedgesinterpsubfcn(x, y(i), X, Y, Z), xmin, xmax, tol);
%             end
            if mode
                lineints(i) = trapz(xi', maxedgesinterpsubfcn(xi', y(i), X, Y, Z));              
            else
                lineints(i) = trapz(xi', zeroedgesinterpsubfcn(xi', y(i), X, Y, Z));
            end
            
            
        end
        
        % Fit a periodic slm object to the periodic line integral data
        ylineintslm = slmengine(y, lineints, ...
            'knots', round(nlineintpnts/3), ...
            'endcon', 'periodic');
        
        interval = [ymin, ymax];
        
    else
        
        if isslm(X) && isscalar(Y) && isscalar(Z)
            ylineintslm = X;
            interval = [Y, Z];
        else
            error('integratehalfperiod2ddata with three inputs expects an slm and ')
        end
        
    end
    
    % Integrate the slm over the interval
    int = periodicslmintegral(ylineintslm, interval);
    
%     if mode
%         int = dblquad(@(x, y) maxedgesinterpsubfcn(x, y, X, Y, Z), xmin, xmax, ymin, ymax, tol, @quadl);
%     else
%         int = dblquad(@(x, y) zeroedgesinterpsubfcn(x, y, X, Y, Z), xmin, xmax, ymin, ymax, tol, @quadl);
%     end
    
end    

function integrand = maxedgesinterpsubfcn(x, y, X, Y, Z)
% here we have a periodic function which is at a maximum at the
% edge of each half period, e.g. cos

    if isoctave
        method = 'cubic';
    else
        method = '*cubic';
    end

    % Test if in odd pole or not
    if round(rem(y-rem(y, 1), 2)) == 0
        % function point is odd, i.e. pole 1, 3 5 etc
        integrand = interp2(X, Y, Z, x(:), repmat(abs(rem(y, 1)), length(x), 1), method);
    else
        % function is even
        integrand = interp2(X, Y, Z, x(:), repmat(1-abs(rem(y, 1)), length(x), 1), method);
    end

end

function integrand = zeroedgesinterpsubfcn(x, y, X, Y, Z)
% here we have a periodic function which crosses the x-axis at the
% edge of each half period, e.g. sin

    if isoctave
        method = 'cubic';
    else
        method = '*cubic';
    end
    
    if round(rem(y-rem(y, 1), 2)) == 0
        % function point is odd, i.e. pole 1, 3 5 etc
        integrand = sign(y) * interp2(X, Y, Z, x(:), repmat(abs(rem(y, 1)), length(x), 1), method);
    else
        % function is even
        integrand = -sign(y) * interp2(X, Y, Z, x(:), repmat(1-abs(rem(y, 1)), length(x), 1), method);
    end

end

% function int = periodicslmintegral(slm, interval)
% 
%     fullperiods = floor(((interval(2) - interval(1)) - (slm.x(end) - slm.x(1)))/(slm.x(end) - slm.x(1)));
%     
%     int = fullperiods * slmpar(slm,'integral');
%         
%     interval = slm.x(1)+mod(interval-slm.x(1), slm.x(end) - slm.x(1));
%     
%     if interval(2) > interval(1)
%     
%         int = int + slmpar(slm,'integral',interval);
%     
%     elseif interval(1) > interval(2)
%         
%         int = int + slmpar(slm,'integral',[slm.x(1),interval(2)]) ...
%                   + slmpar(slm,'integral',[interval(1),slm.x(end)]);
%         
%     end
% 
% end
function int = integratehalfperiodypoly(polymodel, xmin, ymin, xmax, ymax, tol, mode)
% integrateA: finds the double integral of a polynomial fitted to a
% periodic function
%
% Input
%
%   polymodel - a 2d polynomial in x and y whose y values range between 0 and 1
%          where the result of the polynomial between these values is one
%          half period of the actual function to which it is fitted
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
%          integrated, if mode evaluates to true, it is assumed the
%          polynomial is at a maximum at the edges of the domain in the y
%          direction, if false, it is assumed to be zero-crossing at the
%          edges. Default is true for a maximum at the edges.

    if nargin < 6
        tol = 1e-7;
    end
    
    if nargin < 7
        mode = 1;
    end
    
    if mode
        int = dblquad(@(x, y) nestedmaxedgespoly(x, y, polymodel), xmin, xmax, ymin, ymax, tol, @quadl);
    else
        int = dblquad(@(x, y) nestedminedgespoly(x, y, polymodel), xmin, xmax, ymin, ymax, tol, @quadl);
    end
    
end    

function integrand = nestedmaxedgespoly(x, y, polymodel)
% here we have a periodic function which is at a maximum at the
% edge of each half period, e.g. cos

    % Test if in odd pole or not
    if round(rem(y-rem(y, 1), 2)) == 0
        % function point is odd, i.e. pole 1, 3 5 etc
        integrand = polyvaln(polymodel, [x', repmat(abs(rem(y, 1)), length(x), 1)])';
    else
        % function is even
        integrand = polyvaln(polymodel, [x', repmat(1-abs(rem(y, 1)), length(x), 1)])';
    end

end

function integrand = nestedminedgespoly(x, y, polymodel)
% here we have a periodic function which crosses the x-axis at the
% edge of each half period, e.g. sin

    if round(rem(y-rem(y, 1), 2)) == 0
        % function point is odd, i.e. pole 1, 3 5 etc
        integrand = sign(y) * polyvaln(polymodel, [x', repmat(abs(rem(y, 1)), length(x), 1)])';
    else
        % function is even
        integrand = -sign(y) * polyvaln(polymodel, [x', repmat(1-abs(rem(y, 1)), length(x), 1)])';
    end

end
    
    

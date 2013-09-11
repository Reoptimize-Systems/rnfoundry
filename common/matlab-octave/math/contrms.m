function R = contrms(x, Y)
% computes the root mean square of a continuous function sampled at
% intervals
%
% Syntax
%
% R = contrms(x, Y)
%
% Input
%
% Y is a vector or matrix of column vectors of data sampled at the points
% in x, i.e. each column of Y is given by some function
%
% y = f(x)
%
% x is a vector of points at which the values in y were sampled, x should
% be monatonically increasing.
%
% Output
%
% R is the rms of Y (if a vector) or each column of Y if a matrix
%
% See Also CONTMEAN


    if isvector(x) && all(diff(x)>0)
        
        range = x(end) - x(1);

        R = sqrt( trapz(x, Y.^2) / range );
    
    else
        error('x should be a monatonically increasing vector of sample points.');
    end

end
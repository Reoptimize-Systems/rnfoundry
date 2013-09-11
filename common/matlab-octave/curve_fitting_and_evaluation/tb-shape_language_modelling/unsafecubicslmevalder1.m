function y = unsafecubicslmevalder1(x, slm)
% evaluates the first derivative of a cubic slm object with no input
% checking for speed
%
% Syntax
%
% y = unsafecubicslmevalder1(x, slm)
%
% Input
%
%  x - a matrix of independent values at which to evaluate the slm object
%
%  slm - an slm object of degree 3 as produced by slmengine.m
%
% Output
%
%  y - the interpolated values of the first derivative of the slm fit at
%    the positions in x
%

% Created by Richard Crozier 2012

    % clip the points first
%     x = min(slm.knots(end),max(x,knots(1)));

    % and use histc to bin the points.
    [junk,xbin] = histc(x(:),slm.knots);

    % any point which falls at the top end is said to
    % be in the last bin.
    xbin(xbin==slm.nk) = slm.nk-1;

    % first derivative for the cubic case
    t = (x(:)-slm.knots(xbin))./slm.dx(xbin);
    t2 = t.^2;
    s = 1-t;
    s2 = (1-t).^2;
    y = -slm.coef(xbin,2).*(-3*s2+2*s) + ...
        slm.coef(xbin+1,2).*(3*t2-2*t) + ...
        (slm.coef(xbin,1).*(-6*s+6*s2) + ...
        slm.coef(xbin+1,1).*(6*t-6*t2))./slm.dx(xbin);
    
    y = reshape(y, size(x));

end
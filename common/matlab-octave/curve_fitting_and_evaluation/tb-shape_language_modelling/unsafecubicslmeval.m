function y = unsafecubicslmeval(x, slm)
% evaluates a cubic slm object with no input checking for speed
%
% Syntax
%
% y = unsafecubicslmeval(x, slm)
%
% Input
%
%  x - a matrix of independent values at which to evaluate the slm object
%
%  slm - an slm object of degree 3 as produced by slmengine.m
%
% Output
%
%  y - the interpolated values at the positions in x
%

% Created by Richard Crozier 2012

    % clip the points first
%     x = min(slm.knots(end),max(x,knots(1)));

    % and use histc to bin the points.
    [junk,xbin] = histc(x(:),slm.knots);

    % any point which falls at the top end is said to
    % be in the last bin.
    xbin(xbin==slm.nk) = slm.nk-1;

    t = (x(:)-slm.knots(xbin))./slm.dx(xbin);
    t2 = t.^2;
    t3 = t.^3;
    s2 = (1-t).^2;
    s3 = (1-t).^3;
    y = (-slm.coef(xbin,2)   .* (s3 - s2) + ...
          slm.coef(xbin+1,2) .* (t3 - t2)) .* slm.dx(xbin) + ...
          slm.coef(xbin,1)   .* (3*s2 - 2*s3) + ...
          slm.coef(xbin+1,1) .* (3*t2 - 2*t3);
    
    y = reshape(y, size(x));

end
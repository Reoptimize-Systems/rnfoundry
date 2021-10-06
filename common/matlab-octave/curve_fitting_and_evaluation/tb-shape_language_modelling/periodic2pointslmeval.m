function y = periodic2pointslmeval(slm, sep, x, offset)
% evaluates a periodic slm, or two slm objects at two positions either side
% of a specified position
%
% Syntax
%
% y = periodic2pointslmeval(slm, sep, x)
% y = periodic2pointslmeval(slm, sep, x, offset)
%
% Description 
%
% periodic2pointslmeval evaluates the supplied slm object, or array of slm
% objects at two points either side of a position 'x' separated by a
% distance 'sep' such that the slm is evaluated at the points x-sep/2 and
% x+sep/2.
%
% If two slm objects are supplied as an array, the first point (x-sep/2) is
% evaluated using the first slm and the second point (x+sep/2) evaluated
% using the second slm object
%
% An optional offset value can also be supplied which is added to all 'x'
% positions to be evaluated
%
% periodic2pointslmeval returns an (n x 2) matrix of values at the postions
% either side of the supplied n values in x.
%

% Copyright Richard Crozier 2012

    if nargin < 4
        offset = 0;
    end
    
    if numel(slm) == 1
        slminds = [1, 1];
    else
        slminds = [1, 2];
    end
    
    if ~samesizeorscalars(x, offset, sep)
       error('SLMTOOLS:periodic2pointslmeval:badinput', ...
             'x, offset, and sep must be a combination of numeric matrices of the same size or scalar values');
    end

    y = [ periodicslmeval(x(:)+offset-sep(:)/2, slm(slminds(1)), 0, false), ...
          periodicslmeval(x(:)+offset+sep(:)/2, slm(slminds(2)), 0, false) ];
         
end
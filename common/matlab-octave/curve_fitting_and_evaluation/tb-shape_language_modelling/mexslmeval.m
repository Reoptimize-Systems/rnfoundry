% mexslmeval: mex evaluation function for John D'Errico's slm tool
%
% Syntax
%
% y = mexslmeval(x, bins, coef, evalmode)
%
% Input
%
%   x - maxtrix of values of x for which the function is to be interpolated
%
%   bins - vector of bins (the slm.bins field)
%
%   coef - matrix of slm coefficients (the slm.coef field)
%
%   evalmode - numeric flag. Specifies what evaluation is to be done.
%     Unlike slmeval, this is not optional and must be specified:
%
%     == 0 --> evaluate the function at each point in x.
%     == 1 --> the first derivative at each point in x.
%     == 2 --> the second derivative at each point in x.
%     == 3 --> the third derivative at each point in x.
%     == -1 --> NOT YET IMPLEMENTED
%
% Output
%
%    y - Evaluated result, the predicted value, i.e., f(x) or f'(x),
%      f''(x), f'''(x). y will have the same shape and size as x. 
%
% See also: slmeval, slmengine
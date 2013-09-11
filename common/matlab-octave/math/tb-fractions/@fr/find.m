function varargout = find(FRAC,varargin)
% fr/FIND: Find indices of nonzero elements.
%
% See also: ind2sub, non-zeros, relop

% Author: Ben Petschel 28/7/09
%
% Version history:
%   28/7/09 - first release

[varargout{1:nargin}] = find(FRAC~=0,varargin{:});

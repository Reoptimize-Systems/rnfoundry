function answer = ismatnotvec(varargin)
% ISMATRIX      matrix object test
% Returns a vector of logicals (1/0) for which input objects are
% matrices, and which are not (empty, single, vector and arrays). 
%
% answer = ismatrix(varargin)
% 
% See also ISVECTOR ISSCALAR 

    for i=1:length(varargin)
        
      varsize = size(varargin{i});
      
      if length(varsize)<=2 && all(varsize>1)	% varargin{i} is a matrix
        answer(i) = true;			
      else
        answer(i) = false;
      end
      
    end

end
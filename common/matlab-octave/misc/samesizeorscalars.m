function [result, scalarinds, nonscalarinds, nonumericinds] = samesizeorscalars(varargin)
% tests if a set of arguments are a combination of numeric matrices of the
% same size in every dimension and scalar values
% 
% Syntax
%
% [result, scalarinds, nonscalarinds, nonumericinds] = samesizeorscalars(A, B, ...)
%
% Description
%
% samesizeorscalars takes any number of arguments and returns true if they
% are a combination of numeric matrices of the same size in every dimension
% and scalar values
%
% 

% Copyright Richard Crozer 2012

    arescalars = cellfun(@isscalar, varargin);
    
    aresamesize = samesize(varargin{~arescalars});
    
    arenumeric = cellfun(@isnumeric, varargin);
    
    result = every([aresamesize, arenumeric]);
    
    % if requested, get the indices of the scalars and nonscalars etc
    if nargout > 1
        
        scalarinds = find(arescalars);

        if nargout > 2
            
            nonscalarinds = find(~arescalars);
        
            if nargout > 3
                nonumericinds = find(~arenumeric);
            end
            
        end
    
    end
    
end
function result = sameSize(varargin)
% samesize: returns true if dimensions of all the inputs are the same
% 
% Syntax
% 
% result = check.sameSize(A, B, ...)
% 
% Description
% 
% A, B etc are arrays or matrices of any type. samesize returns true if
% they are all the same size in every dimension, and have the same number
% of dimensions, or false otherwise
% 
% check.sameSize returns true for zero arguments, or a single argument on
% the basis that the input will definately be the same size as itself in
% these cases.
% 
% Example
% 
% astruct = struct('a',{1,2,3}); % creates a (1 x 3) struct array
% anarray = [1,2,3];
% acellarray = {2+6j, [1,2,3,4,5,6,7,8], 'a string'};
% check.sameSize(astruct, anarray, acellarray)
% 
% >>
%  ans =
%       1
% 
% check.sameSize(astruct, anarray', acellarray)
% 
% >>
%  ans =
%       0
% 
% See also: size, ndims
% 
     
% Created 2018 Richard Crozer

    % return true if less than two arguments are passed in
    if nargin < 2
        result = true;
        return;
    end
    
    % now test the actual sizes of the each input for equality
    dimscells = cellfun(@size, varargin, 'UniformOutput', false);

    % if all have the same number of dimensions, return true if the sizes
    % of all dimensions are the same or false otherwise
    result = all(isequal(dimscells{:}));

end
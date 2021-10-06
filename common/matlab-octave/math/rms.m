function rmsval = rms(x)
% rms: calculates the rms (root mean square) value of a vector
%
% Input:
%
%   x - vector of values for which the rms value is to be calculated
%
% Output:
%
%   rmsval - the rms value of the numbers in the vector

    % norm(x) when x is a vector returns the following formula:
    %
    % sum(abs(x).^2)^(1/2)
    %
    % the formula for rms is:
    %
    % sqrt(sum(x.^2) ./ n)
    %
    % where n is the number of samples, leading to the following statement:
    
    rmsval = norm(x) / sqrt(length(x));

end
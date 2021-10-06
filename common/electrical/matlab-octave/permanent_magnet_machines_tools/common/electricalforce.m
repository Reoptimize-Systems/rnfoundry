function F = electricalforce(I, dlambdadx)
% electricalforce: calculates the force from an on-load electrical machine
% from the current and rate of change of flux linkage with respect to
% displacement at the same instant.

    F = sum(I .* dlambdadx);

end
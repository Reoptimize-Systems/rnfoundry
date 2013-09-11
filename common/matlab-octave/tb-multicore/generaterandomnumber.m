function [randomNr, randomStr] = generaterandomnumber
%GENERATERANDOMNUMBER
%   in very unlikely cases, it might happen that the random states of rand
%   and randn are equal in two Matlab processes calling function
%   SETFILESEMAPHORE. For this reason, the system and cpu time are used to
%   create different random numbers even in this unlikely case.
%
%	This all were not necessary if it were possible to get some sort of a
%	Matlab process ID.

    nOfDigits = 8; % length of random string will be 4*nOfDigits

    randNr    = rand;
    randnNr   = mod(randn+0.5, 1);
    cputimeNr = mod(cputime, 100)/100;
    nowNr     = mod(rem(now,1)*3600*24, 100)/100;

    % random number is used for random pause after conflict
    randomNr = 0.25 * (randNr + randnNr + cputimeNr + nowNr);

    % random string is used for the semaphore file name
    if nargout > 1
        ee = 10^nOfDigits;
        randomStr = sprintf('%.0f%.0f%.0f%.0f', ...
            ee * randNr,    ...
            ee * randnNr,   ...
            ee * cputimeNr, ...
            ee * nowNr      ...
            );
    end

end
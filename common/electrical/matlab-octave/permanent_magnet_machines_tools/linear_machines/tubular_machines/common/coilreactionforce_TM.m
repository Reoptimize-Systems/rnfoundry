function shearForce = coilreactionforce_TM(RoVRm, Rm, g, N, avRadBvals, I, div)
% Determines the shear force ue to the interation of the current carrying
% coils and the radial magnetic field in a tubular machine
%
% Input:
%
%   RoVRm - Ro \ Rm ratio
%
%   Rm - the translator radius in m
%
%   g - the airgap in m
%
%   N - total number of turns in the winding
%
%   avRadBvals - (n x p) matrix of the average radial flux density in the
%                coil sections. Each column in avRadBvals corresponds to a
%                different coil position.
%
%   I - scalar or (1 x p) vetor of values of the current in the coil at
%       each of the positions for which the average radial flux density
%       values are given in avRadBvals. If a scalar, the same value is
%       assumed for all positions
%
%   div - the number of segments the coil is split into for the purposes of
%   analysing the net force on a translator displaced from centre
%
% Output:
%
%   shearForce - Row vector of values of the shear force on the coil at
%                each position corresponding to a column in avRadBvals
%                

    Stot = size(avRadBvals, 1);
    
    % Determine the total length of wire in the winding
    length = wirelength_TM(RoVRm, N, Rm, g);
    
    % Get length of one turn
    length = length / N; 
    
    % get length of wire to be considered for force calculation
    delLength = length / div;
    
    % For each of the sets of average radial flux densities, each one
    % corresponding to a different coil position
    shearForce = sum(avRadBvals,1) .* I .* delLength .* (N/Stot);

end
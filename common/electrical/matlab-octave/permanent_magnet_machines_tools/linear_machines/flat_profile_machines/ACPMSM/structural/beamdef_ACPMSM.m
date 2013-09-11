function Def = beamdef_ACPMSM(IVars, totalLength, aL, x, FVars, E, IMethod, forceRatio)
% function: beamdef_ACPMSM
% 
% Calculates the deflection of the parts of the air-cored pmsm with
% some distributed loading using
%
% Input: 
%
%   IVars - (1 x 4) matrix of values for calculating the second moment of
%           area of the stator and translator sections, the second row of
%           values is used for both stator sections if a double sided
%           machine is being investigated.
%
%           IVars(1,1): b, width of the I-beam flanges
%           IVars(1,2): t, the thickness of the flanges
%           IVars(1,3): tw, the thickness of the I-Beam vertical part
%           IVars(1,4): d, height of the I-Beam vertical part
%
%   totalLength - total length of the beam
%
%   aL - Either: 1. If scalar, the 'active' length of translator, i.e the length of
%                beam undergoing the loading.
%
%                2. If FVars's length is greater than 2, aL must be an array
%                (of size one more than FVars) containing the positions at
%                which the forces described in FVars are applied. In this
%                case FVars values are the total applied between each adjacent
%                value of aL, i.e. the sum of the forces between aL(n) and
%                aL(n+1) would be in FVars(n).
%
%   x - vector of position values at which the deflection is to be
%       calculated 
% 
%   FVars - Either: 1. If row vector, a uniform force acting on the active
%                   length of each part of the machine.
%
%                   2. If (n x 2) matrix, a linearly distributed force along
%                   the active length for each part.
%
%                   2. If (n x p) matrix, each value in a row of FVars 
%                   are the total force applied between each adjacent value
%                   of aL, i.e. the sum of the forces between aL(n) and
%                   aL(n+1) for each part of the machine.
%
%   E - Scalar containing Young's modulus of the I-beam material.
%
% Output:
%
%   Def - values of the deflections at each x position for each part of the
%         machine in question 
%
    
    if nargin < 6
       % If no modulus of elasticity supplied, use 200 GPa
       % steel at room temp is 207 GPa 
       E = 207e9;
    end
    
    % Here we define both deflections and forces upward to be positive, but forces
    % downwards are defined as positive in the equations, therefore we must
    % flip the sign of the forces in the input
    FVars = -FVars;

    for i = 1:size(FVars,1)
        
        switch size(FVars,2)
            
            case 0
                error('No forces supplied')

            case 1
                % uniform force
                a2 = totalLength - ((totalLength - aL)/2);
                a1 = (totalLength - aL)/2;
                Yvars(1,:,i) = [FVars(i,1) FVars(i,1) totalLength a1];
                if a2 == totalLength
                    % Supply dummy Yvars variables as a value of a = l
                    % (length) results in an attempt to divide by zero.
                    Yvars(2,:,i) = [0 0 totalLength 0];
                else
                    Yvars(2,:,i) = [-FVars(i,1) -FVars(i,1) totalLength a2];
                end
            case 2
                % linearly variable force
                a2 = totalLength - ((totalLength - aL)/2);
                a1 = (totalLength - aL)/2;
                wm = FVars(i,2);
                wa = FVars(i,1);
                Yvars(:,:,i) = CalculateYvarsDistribF(a1, a2, totalLength, wa, wm);

            otherwise
                
                if (size(aL,2)-1) == size(FVars,2)
                    % multiple distributed forces, where the forces in
                    % FVars are the TOTAL force applied between each adjacent
                    % value of aL, describing the location of the forces along
                    % the beam, i.e. the sum of the forces between aL(n) and
                    % aL(n+1).
                    Yvars(:,:,i) = CalculateYvarsDistribF(aL(1:(end-1))', aL(2:end)', totalLength, FVars(i,:)');
                else
                    error('Too many forces supplied') 
                end 
        end
        
    end
    
    beamMethod = '3.2d';
    
    % Calculate deflections in magnet supports using I-Beam method
    % (1.6)
    if nargin < 7
       IMethod = '1.6' ;
    end
    
    if nargin < 8
        forceRatio = 1;
    end
    
    Def(1,:) = BeamDeflectionSuperY1([Yvars(:,1:2,1) .* forceRatio, Yvars(:,3:end,1)], ...
                                        IVars(1,:), E(1), x, IMethod, beamMethod);

    Def(2,:) = BeamDeflectionSuperY1([Yvars(:,1:2,2) .* forceRatio, Yvars(:,3:end,2)], ...
                                        IVars(1,:), E(1), x, IMethod, beamMethod);
  
end
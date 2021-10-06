function Def = beamdef_PMSM(IVars, totalLength, aL, x, FVars, E, IMethod, forceRatio, beamMethod)
% Calculates the deflection of the parts of the pmsm with some distributed
% loading using
%
% Input: 
%
%   IVars - (1 x 2) cell array of values, each of which contains a matrix
%           of values for calculating the second moment of area of the
%           stator and translator sections repectively, the the values in
%           the second cell are used for both stator sections if a double
%           sided machine is being investigated.
%
%           The contents of the first cell of IVars is determined by the
%           type of beam to be used to support the structur, as specified
%           in the IMethod field. 
%
%           If an I-Beam is to be used it will contain the following values
%
%           IVars{1}(1,1): b, width of the I-beam flanges
%           IVars{1}(1,2): t, the thickness of the flanges
%           IVars{1}(1,3): tw, the thickness of the I-Beam vertical part
%           IVars{1}(1,4): d, height of the I-Beam vertical part
%
%           If an alternative section is to be used, it will contain
%           appropriate values for that beam.
%
%           The second cell of IVars
%
%           IVars{2}(1,1): dbi, half thickness of central section (back 
%                          iron thickness)
%           IVars{2}(1,2): Taup, the pole width
%           IVars{2}(1,3): hs, the tooth height
%           IVars{2}(1,4): bt, the tooth width
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
%   E - (1 x 2) Vector containing Young's modulus of the beam materials.
%       The first value is Young Modulus of the I-beam support on the
%       stator(s) the second is for the translator.
%
%   forceRatio - If the translator support beam has a different surface
%                area to the armature 'beam', we must adjust the forces on
%                the beams to account for this and calculate the correct
%                deflections (the support beams will be spaced with no
%                regard for pole size). forceRatio holds the ratio of the
%                support beam width to the pole width. If not supplied, the
%                default value of 1 is used.
%
% Output:
%
%   Def - values of the deflections at each x position for each part of the
%         machine in question 
%
    
    if nargin < 6
       % If no modulus of elasticity supplied, use 200 GPa
       % steel at room temp is 207 GPa 
       E = [200e9 151e9]; 
       
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
    
    
    if nargin < 7
        IMethod = '1.6';
    end
    
    if nargin < 8
        forceRatio = 1;
    end
    
    if nargin < 9
        beamMethod = '3.2d';
    end
    
    
    if size(FVars,1) == 2
        % First calculate deflection in stator support using I-Beam method
        % (1.6). We multiply the forces by the force ratio calculated
        % earlier to get the correct force per unit length
        Def(1,:) = BeamDeflectionSuperY1([Yvars(:,1:2,1) .* forceRatio, Yvars(:,3:end,1)], IVars{1}(1,:), E(1), x, IMethod, beamMethod);
        % Next calculate deflection in translator, passing in function name
        % as it's not a standard formula
        Def(2,:) = BeamDeflectionSuperY1(Yvars(:,:,2), IVars{2}(1,:), E(2), x, 'TranslatorMomentOfArea_SS_PMSM', beamMethod);
        
    elseif size(FVars,1) == 3
        % First calculate deflection in upper stator support using I-Beam method
        % (1.6). We multiply the forces by the force ratio calculated
        % earlier to get the correct force per unit length
        Def(1,:) = BeamDeflectionSuperY1([Yvars(:,1:2,1) .* forceRatio, Yvars(:,3:end,1)], IVars{1}(1,:), E(1), x, IMethod, beamMethod);
        % Next calculate deflection in translator, passing in function name
        % as it's not a standard formula
        Def(2,:) = BeamDeflectionSuperY1(Yvars(:,:,2), IVars{2}(1,:), E(2), x, 'TranslatorMomentOfArea_DS_PMSM', beamMethod);
        % Next calculate deflection in lower stator support using I-Beam method
        % (1.6). Again we multiply the forces by the force ratio.
        Def(3,:) = BeamDeflectionSuperY1([Yvars(:,1:2,3) .* forceRatio, Yvars(:,3:end,3)], IVars{1}(1,:), E(1), x, IMethod, beamMethod);
        
    end
       
end
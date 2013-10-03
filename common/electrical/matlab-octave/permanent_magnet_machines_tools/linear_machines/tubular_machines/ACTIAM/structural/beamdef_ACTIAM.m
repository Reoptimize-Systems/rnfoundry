function Def = beamdef_ACTIAM(IVars, supportLengths, Lf, La, P, fFVars, aFVars, E, x)
% function: beamdef_ACTIAM
% 
% Calculates the deflection of the parts of the air-cored tubular machine
% with an outer sheath with some distributed loading using euler beam
% theory
%
% Input: 
%
%   IVars - (3 x 2) matrix of values for calculating the second moment of
%           area of the stator and translator sections, the second row of
%           values is used for both stator sections if a double sided
%           machine is being investigated.
%
%           IVars(1,1): Rso, outer radius of shaft support
%           IVars(1,2): Rsi, inner radius of shaft support
%
%           IVars(2,1): Ro, outer coil radius
%           IVars(2,2): Ri, inner coil radius
%
%           IVars(3,1): Ra, outer sheath radius
%           IVars(3,2): Ro, inner sheath radius
%
%   supportLengths - (1 X 2) vector of lengths of the supports at either
%                    end of the active part of the field/translator. I.e.
%                    the first value is the distance from the first bearing
%                    to the start of the active part of the translator,
%                    while the second is the distance from the end of the
%                    active part to the second bearing.
%
%   Lf - Either: 1. If scalar, the 'active' length of translator, i.e the length of
%                beam undergoing the loading.
%
%                2. If FVars's length is greater than 1, Lf must be an array
%                (of size one more than FVars) containing the positions at
%                which the forces described in FVars are applied. In this
%                case FVars values are the total applied between each adjacent
%                value of Lf, i.e. the sum of the forces between Lf(n) and
%                Lf(n+1) would be in FVars(n).
%
%   P - scalar, the total axial loading acting on the machine
%
%   fFVars - Either: 1. If row vector, a uniform force acting on the active
%                   length of each part of the machine.
%
%                   2. If (n x 2) matrix, a linearly distributed force along
%                   the active length for each part.
%
%                   2. If (n x p) matrix, each value in a row of FVars 
%                   are the total force applied between each adjacent value
%                   of Lf, i.e. the sum of the forces between Lf(n) and
%                   Lf(n+1) for each part of the machine.
%
%   E - (1 x 3) Vector containing Young's modulus of the beam materials.
%       The first value is Young Modulus of the shaft, the second the the
%       young's modulus of the coil and the third the young's modulus of
%       the outer sheath
%
%   x - vector of position values at which the deflection is to be
%       calculated 
%
% Output:
%
%   Def - values of the deflections at each x position for each part of the
%         machine in question 
%
    
    if nargin < 6
       % If no modulus of elasticity supplied, use 200 GPa
       % steel at room temp is 207 GPa 
       E = [200e9 200e9]; 
       
    end
    
    
    % Here we define both deflections and forces upward to be positive, but forces
    % downwards are defined as positive in the equations, therefore we must
    % flip the sign of the forces in the input
    fFVars = -fFVars;
    aFVars = -aFVars;
    
    switch size(fFVars,2)

        case 0
            error('No forces supplied')

        case 1
            % uniform force
            fYvars = [fFVars(1,1) fFVars(1,1) Lf(end) 0];
        case 2
            % linearly variable force
            wm = fFVars(1,2);
            wa = fFVars(1,1);
            fYvars = CalculateYvarsDistribF(0, Lf, Lf, wa, wm);

        otherwise

            if (size(Lf,2)-1) == size(fFVars,2)
                % multiple distributed forces, where the forces in
                % FVars are the TOTAL force applied between each adjacent
                % value of Lf, describing the location of the forces along
                % the beam, i.e. the sum of the forces between Lf(n) and
                % Lf(n+1).
                fYvars = CalculateYvarsDistribF(Lf(1:(end-1))', Lf(2:end)', totalLength, fFVars(1,:)');
            else
                error('Too many forces supplied') 
            end 
    end
    
    switch size(aFVars,2)

        case 0
            error('No forces supplied')

        case 1
            % uniform force
            aYvars(1,:) = [aFVars(1,1) aFVars(1,1) La supportLengths(1)];
            aYvars(2,:) = [-aFVars(1,1) -aFVars(1,1) La (La - supportLengths(2))];
        case 2
            % linearly variable force
            wm = aFVars(1,2);
            wa = aFVars(1,1);
            aYvars = CalculateYvarsDistribF(supportLengths(1), La - supportLengths(2), La, wa, wm);

        otherwise
            
            if (size(La,2)-1) == size(aFVars,2)
                % multiple distributed forces, where the forces in
                % aFVars are the TOTAL force applied between each adjacent
                % value of La, describing the location of the forces along
                % the armature, i.e. the sum of the forces between La(n) and
                % La(n+1).
                aYvars = CalculateYvarsDistribF(La(1:(end-1))', La(2:end)', totalLength, aFVars(1,:)');
            else
                error('Too many forces supplied') 
            end 
    end
    
    % Add the axial loads to the Yvars matrices
    fYvars = [repmat(P, size(fYvars,1), 1) fYvars];
    aYvars = [repmat(P, size(aYvars,1), 1) aYvars];
    
    beamMethod = '10.2f'; % Left end guided, right end simply supported, axial load
    
    % The forces on the field will be in fFvars, 
    Def = BeamDeflectionSuper(fYvars, IVars(1,:), E(1), x, '1.15', beamMethod);
    
    % The sheath and coils will act as one member but are of different
    % materials, we must therefore create an equivalent beam of one
    % material by changing the moment of area of either the sheath or
    % coils. We will do this by making a new beam of the appropriate
    % cross-sectional area.
    Icoil = MomentOfArea(IVars(2,:), '1.15');
    
    % Calculate the I we require if the coil were made of the same material
    % as the sheath
    Icoil = Icoil * (E(2) / E(3));
    
    % Determine the new IVars we require to achieve this I for the combined
    % sheath and coil beam
    modifiedIVars = combinedsheatandcoilI_TM(Icoil, IVars(2:3,:));

    % Next calculate deflection in the armature,  using appropriate
    % method as stored in armatureIMethod
    Def(2,:) = BeamDeflectionSuper(aYvars, modifiedIVars, E(4), x, '1.15', beamMethod);    
         
end

%% Old Code (cantilever supports)
%   IVars - (3 x 2) matrix of values for calculating the second moment of
%           area of the stator and translator sections, the second row of
%           values is used for both stator sections if a double sided
%           machine is being investigated.
%
%           IVars(1,1): Rso, outer radius of shaft support
%           IVars(1,2): Rsi, inner radius of shaft support
% 
%           IVars(2,1): Rm, outer radius of field
%           IVars(2,2): Rb, radius of bore hole in field for tension cable
%
%           IVars(3,1): Ro, outer coil radius
%           IVars(3,2): Ri, inner coil radius
%
%           IVars(4,1): Ra, outer sheath radius
%           IVars(4,2): Ro, inner sheath radius
%

%     % we will build the YVars
%     % Matrix for the field supports from the sum of these
%     sFVars = sum(fFVars(1,:),2)/2;
%     
%     sYvars = [sFVars supportLengths(1) 0; sFVars supportLengths(2) 0];
%     
%     beamMethod = '3.1f'; % Use cantilever
%     % Calculate deflection and slope in field support clamped to end of active part
%     % at it's end, we will use formula for thick cylinder
%     supportDef = BeamDeflection(IVars(1,1:2), sYvars(1,:), E(1), supportLengths(1), '1.15', beamMethod);
%     
%     supportSlope = BeamSlope(IVars(1,1:2), sYvars(1,:), E(1), supportLengths(1), '1.15', beamMethod);
%     
%     supportDef(2) = BeamDeflection(IVars(1,1:2), sYvars(2,:), E(1), supportLengths(2), '1.15', beamMethod);
%     
%     supportSlope(2) = BeamSlope(IVars(1,1:2), sYvars(2,:), E(1), supportLengths(2), '1.15', beamMethod);
%     
%     % Find the additional deflection  in the active part due to the deflection of the supports
%     initY = supportDef(1) - ((x(x>supportLengths(1)&x<(La(end)-supportLengths(2)))-supportLengths(1)) .* sin(asin((supportDef(1)-supportDef(2))/Lf(end))));
%     
%     % Calculate displacement at x positions in cantilevers and insert into
%     % deflection matrix at the appropriate points
%     Def = BeamDeflection(IVars(1,1:2), sYvars(1,:), E(1), fliplr(x(x<=supportLengths(1))), '1.15', beamMethod);
%     
%     Def(1,size(Def,2)+1:size(x,2)) = 0;
%     
%     Def(1,size(x,2)-size(x(x<supportLengths(2)),2)+1:size(x,2)) = BeamDeflection(IVars(1,1:2), sYvars(2,:), E(1), fliplr(La(end) - x(x>=(La(end)-supportLengths(2)))), '1.15', beamMethod);
% 
%     % Calculate deflection in the field
%     beamMethod = '3.2d'; % Use fixed ends
%     
%     % Calculate the relative angular displacement of the beam from the two
%     % supportSlope values
%     thetaYVars = [ (supportSlope(1)-sin(asin((supportDef(1)-supportDef(2))/Lf(end)))) Lf(end) 0;...
%         (supportSlope(2)+sin(asin((supportDef(1)-supportDef(2))/Lf(end)))) Lf(end) Lf(end)];
%     
%     Def(1,size(x(x<=supportLengths(1)),2)+1:size(x(x<(La(end)-supportLengths(2))),2)) = -BeamDeflectionSuper(fYvars, IVars(2,1:2), E(2), x(x>supportLengths(1)&x<(La(end)-supportLengths(2)))-supportLengths(1), '1.15', beamMethod)...
%         - BeamDeflectionSuper(thetaYVars, IVars(2,1:2), E(2), x(x>supportLengths(1)&x<(La(end)-supportLengths(2)))-supportLengths(1), '1.15', '3.4d')...
%         - initY;
%
%     % The sheath and coils will act as one member but are of different
%     % materials, we must therefore create an equivalent beam of one
%     % material by changing the moment of area of either the sheath or
%     % coils. We will do this by making a new beam of the appropriate
%     % cross-sectional area.
%     Icoil = MomentOfArea(IVars(3,:), '1.15');
%     
%     % Calculate the I we require if the coil were made of the same material
%     % as the sheath
%     Icoil = Icoil * (E(3) / E(4));
%     
%     % Determine the new IVars we require to achieve this I for the combined
%     % sheath and coil beam
%     modifiedIVars = combinedsheatandcoilI_TM(Icoil, IVars(3:4,:));
% 
%     % Next calculate deflection in the armature,  using appropriate
%     % method as stored in armatureIMethod
%     Def(2,:) = BeamDeflectionSuper(aYvars, modifiedIVars, E(4), x, '1.15', beamMethod);    
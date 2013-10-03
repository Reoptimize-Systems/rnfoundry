function varargout = beamdef_TM(IVars, supportLengths, totalLength, aL, P, fFVars, aFVars, x, E, suppWperm, beamMethod)
% function: beamdef_TM
% 
% Calculates the deflection of the parts of the tubular machine
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
%   supportLengths - (2 X 2) matrix of lengths of the supports at either
%                    end of the active part of the field/translator. I.e.
%                    the first value is the distance from the first bearing
%                    to the start of the active part of the translator,
%                    while the second is the distance from the end of the
%                    active part to the second bearing.
%
%   totalLength - (1 X 2) Vector. The total length of the field and
%                 armature support respectively.
%
%   aL - Either: 1. If scalar, the 'active' length of armature, i.e the length of
%                the part of both beams undergoing the loading.
%
%                2. If FVars's length is greater than 1, aL must be an array
%                (of size one more than FVars) containing the positions at
%                which the forces described in FVars are applied. In this
%                case FVars values are the total applied between each adjacent
%                value of Lf, i.e. the sum of the forces between aL(n) and
%                aL(n+1) would be in FVars(n).
%
%   P - scalar, the total axial loading acting on the machine, the total of
%       the shear forces
%
%   fFVars - Either: 1. If row vector, a uniform force acting on the active
%                   length of each part of the machine.
%
%                   2. If (n x 2) matrix, a linearly distributed force along
%                   the active length for each part.
%
%                   2. If (n x p) matrix, each value in a row of FVars 
%                   are the total force applied between each adjacent value
%                   of aL, i.e. the sum of the forces between aL(n) and
%                   aL(n+1) for each part of the machine. This does not
%                   include the forces on the supports.
%
%   aFVars - Either: 1. If row vector, a uniform force acting on the active
%                   length of each part of the machine.
%
%                   2. If (n x 2) matrix, a linearly distributed force along
%                   the active length for each part.
%
%                   2. If (n x p) matrix, each value in a row of FVars 
%                   are the total force applied between each adjacent value
%                   of aL, i.e. the sum of the forces between aL(n) and
%                   aL(n+1) for each part of the machine. This does not
%                   include the forces on the supports.
%
%   x - vector of position values at which the deflection is to be
%       calculated 
%
%   E - (1 x 3) Vector containing Young's modulus of the beam materials.
%       The first value is Young Modulus of the shaft, the second the the
%       young's modulus of the coil and the third the young's modulus of
%       the outer sheath [shaft coil sheath]
%
%   suppWperm - optional (1 x 2) vector of weight per m on the supports of
%               the field and armature respectively for including their own
%               weight in the analysis, can be omitted for finding
%               incremental forces
%
% Output:
%
%   Def - values of the deflections at each x position for each part of the
%         machine in question 
%
    
    if nargin < 9
       % If no modulus of elasticity supplied, use 200 GPa
       % steel at room temp is 207 GPa 
       E = [200e9 100 200e9]; 
       
    end
     
    switch size(fFVars,2)

        case 0
            error('No forces supplied')

        case 1
            % Uniform Force
            % First determine the position of the end of the active part
            a2 = totalLength(1) - supportLengths(1,2);
            % Next determine the position of the start of the active part
            a1 = supportLengths(1,1);
            % Create a vector describing a uniform loading from the end of
            % the beam to the point a1 at the start of the active part
            fYvars(1,:) = [fFVars(1,1) fFVars(1,1) totalLength(1) a1];
            if a2 == totalLength(1)
                % If the active part runs right to the end of the support,
                % supply dummy Yvars variables as a value of a = l
                % (length) results in an attempt to divide by zero.
                fYvars(2,:) = [0 0 totalLength(1) 0];
            else
                % Otherwise, create a vector describing an upward uniform
                % loading from the end of the beam to the point a2 at the
                % start of the active part to cancel out the forces at this
                % part
                fYvars(2,:) = [-fFVars(1,1) -fFVars(1,1) totalLength(1) a2];
            end
            
            
        case 2
            % linearly variable force
            a2 = totalLength - supportLengths(1,2);
            a1 = supportLengths(1,1);
            wm = fFVars(1,2);
            wa = fFVars(1,1);
            fYvars(:,:) = CalculateYvarsDistribF(a1, a2, totalLength(1), wa, wm);

        otherwise

            if (size(aL,2)-1) == size(fFVars,2)
                % multiple distributed forces, where the forces in
                % FVars are the TOTAL force applied between each adjacent
                % value of aL, describing the location of the forces along
                % the beam, i.e. the sum of the forces between aL(n) and
                % aL(n+1).
                if supportLengths(1,1) > 0 && supportLengths(1,2) > 0
                    tempaL = [aL + supportLengths(1,1), aL(end)+supportLengths(1,1)+supportLengths(1,2)];
                    a1 = [0 tempaL(1:(end-1))]';
                    a2 = tempaL';
                    wa = [0 fFVars 0]';
                elseif supportLengths(1,1) > 0
                    tempaL = aL + supportLengths(1,1);
                    a1 = [0 tempaL(1:(end-1))]';
                    a2 = tempaL';
                    wa = [0 fFVars]';
                elseif supportLengths(1,2) > 0
                    tempaL = [aL + supportLengths(1,1), aL(end)+supportLengths(1,1)+supportLengths(1,2)];
                    a1 = tempaL(1:(end-1))';
                    a2 = tempaL(2:end)';
                    wa = [fFVars 0]';
                else
                    a1 = aL(1:(end-1))';
                    a2 = aL(2:end)';
                    wa = fFVars';
                end
                
                fYvars = CalculateYvarsDistribF(a1, a2, totalLength(1), wa);

            else
                error('Too many forces supplied') 
            end 
    end
    
    if nargin == 10
        if size(suppWperm,2) ~= 2
            error('If supplying support weights, you must supply them for both the field and armature')
        end
        % First determine the position of the end of the active part
        a2 = totalLength(1) - supportLengths(1,2);
        % Next determine the position of the start of the active part
        a1 = supportLengths(1,1);
        % If the support weights are supplied create a loading case
        % to include their own weight
        fYvars(end+1,:) = [suppWperm(1,1) suppWperm(1,1) totalLength(1) 0];
        fYvars(end+1,:) = [-suppWperm(1,1) -suppWperm(1,1) totalLength(1) a1];

        if a2 == totalLength(1)
            % If the active part runs right to the end of the support,
            % supply dummy Yvars variables as a value of a = l
            % (length) results in an attempt to divide by zero.
            fYvars(end+1,:) = [0 0 totalLength(1) 0];
        else
            % Otherwise, create a vector describing an upward uniform
            % loading from the end of the beam to the point a2 at the
            % start of the active part to cancel out the forces at this
            % part
            fYvars(end+1,:) = [suppWperm(1,1) suppWperm(1,1) totalLength(1) a2];
        end
    end
    
    switch size(aFVars,2)

        case 0
            error('No forces supplied')

        case 1
            % uniform force
            a2 = totalLength(2) - supportLengths(2,2);
            a1 = supportLengths(2,1);
            aYvars(1,:) = [aFVars(1,1) aFVars(1,1) totalLength(2) a1];
            if a2 == totalLength(2)
                % Supply dummy Yvars variables as a value of a = l
                % (length) results in an attempt to divide by zero.
                aYvars(2,:) = [0 0 totalLength(2) 0];
            else
                aYvars(2,:) = [-aFVars(1,1) -aFVars(1,1) totalLength(2) a2];
            end
        case 2
            % linearly variable force
            a2 = totalLength(2) - supportLengths(2,2);
            a1 = supportLengths(2,1);
            wm = aFVars(1,2);
            wa = aFVars(1,1);
            aYvars = CalculateYvarsDistribF(a1, a2, totalLength(2), wa, wm);

        otherwise

            if (size(aL,2)-1) == size(aFVars,2)
                % multiple distributed forces, where the forces in
                % FVars are the TOTAL force applied between each adjacent
                % value of aL, describing the location of the forces along
                % the beam, i.e. the sum of the forces between aL(n) and
                % aL(n+1).
                if supportLengths(2,1) > 0 && supportLengths(2,2) > 0
                    tempaL = [aL + supportLengths(2,1), aL(end)+supportLengths(2,1)+supportLengths(2,2)];
                    a1 = [0 tempaL(1:(end-1))]';
                    a2 = tempaL';
                    wa = [0 aFVars 0]';
%                     aYvars = CalculateYvarsDistribF([0 aL(1:(end-1))]', aL', totalLength(2), [0 aFVars 0]');
                elseif supportLengths(2,1) > 0
                    tempaL = aL + supportLengths(2,1);
                    a1 = [0 tempaL(1:(end-1))]';
                    a2 = tempaL';
                    wa = [0 aFVars]';
%                     aYvars = CalculateYvarsDistribF([0 aL(1:(end-1))]', aL', totalLength(2), [0 aFVars]');
                elseif supportLengths(2,2) > 0
                    tempaL = [aL + supportLengths(2,1), aL(end)+supportLengths(2,1)+supportLengths(2,2)];
                    a1 = tempaL(1:(end-1))';
                    a2 = tempaL(2:end)';
                    wa = [aFVars 0]';
                else
                    a1 = aL(1:(end-1))';
                    a2 = aL(2:end)';
                    wa = aFVars';
                end
                
                aYvars = CalculateYvarsDistribF(a1, a2, totalLength(2), wa);
                
            else
                error('Too many forces supplied') 
            end 
    end
    
    if nargin == 10
        if size(suppWperm,2) ~= 2
            error('If supplying support weights, you must supply them for both the field and armature')
        end
        % First determine the position of the end of the active part
        a2 = totalLength(1) - supportLengths(2,2);
        % Next determine the position of the start of the active part
        a1 = supportLengths(2,1);
        % If the support weights are supplied create a loading case
        % to include their own weight
        aYvars(end+1,:) = [suppWperm(2) suppWperm(2) totalLength(1) 0];
        aYvars(end+1,:) = [-suppWperm(2) -suppWperm(2) totalLength(1) a1];

        if a2 == totalLength(1)
            % If the active part runs right to the end of the support,
            % supply dummy Yvars variables as a value of a = l
            % (length) results in an attempt to divide by zero.
            aYvars(end+1,:) = [0 0 totalLength(1) 0];
        else
            % Otherwise, create a vector describing an upward uniform
            % loading from the end of the beam to the point a2 at the
            % start of the active part to cancel out the forces at this
            % part
            aYvars(end+1,:) = [suppWperm(2) suppWperm(2) totalLength(1) a2];
        end
    end
    
    % Add the axial loads to the Yvars matrices, this must be added to
    % every case
    fYvars = [repmat(P, size(fYvars,1), 1) fYvars];
    aYvars = [repmat(P, size(aYvars,1), 1) aYvars];
    
    if nargin < 11
        beamMethod = '10.2e'; % Left end simply supported, right end simply supported, axial load
    end
    
    % The forces on the field will be in fFvars, and the x coordinates will
    % agree
    Def = BeamDeflectionSuperY1(fYvars, IVars(1,:), E(1), x, '1.15', beamMethod);
    
    % We can also determine the moments in the beam in order to determine
    % whether we are in the linear region of the material's stress-strain
    % curve
    Mom = BeamMomentSuperY1(fYvars, IVars(1,:), E(1), x, '1.15', beamMethod);
    
    % The sheath and coils will act as one member but are of different
    % materials, we must therefore create an equivalent beam of one
    % material by changing the moment of area of either the sheath or
    % coils. We will do this by making a new beam of the appropriate
    % cross-sectional area.
    Icoil = MomentOfInertiaY1(IVars(2,:), '1.15');
    
    % If the coil section or support is shorter than the shaft length, we
    % must ensure the x coordinates are correct and that the deflections
    % are mapped to the correct positions
    if totalLength(2) < totalLength(1)
       x = x - ((totalLength(1) - totalLength(2))/2);
       defshift = size(x(x<0),2)+1;
       x = x(x>=0 & x<=totalLength(2));
    else
        defshift = 1;
    end
    
    Def(2,:) = zeros(1,size(Def,2));
    
    if size(E,2) == 3
        % Calculate the I we require if the coil were made of the same material
        % as the sheath
        Icoil = Icoil * (E(2) / E(3));

        % Determine the new IVars we require to achieve this I for the combined
        % sheath and coil beam
        modifiedIVars = combinedsheatandcoilI_TM(Icoil, IVars(2:3,:));

        % Next calculate deflection in the armature, using appropriate
        % method as stored in armatureIMethod
        Def(2,defshift:(defshift+size(x,2)-1)) = BeamDeflectionSuperY1(aYvars, modifiedIVars, E(3), x, '1.15', beamMethod); 
        Mom(2,defshift:(defshift+size(x,2)-1)) = BeamMomentSuperY1(aYvars, modifiedIVars, E(3), x, '1.15', beamMethod);
    else
        % Next calculate deflection in the armature,  using appropriate
        % method as stored in armatureIMethod
        Def(2,defshift:(defshift+size(x,2)-1)) = BeamDeflectionSuperY1(aYvars, IVars(2,:), E(2), x, '1.15', beamMethod);
        Mom(2,defshift:(defshift+size(x,2)-1)) = BeamMomentSuperY1(aYvars, IVars(2,:), E(2), x, '1.15', beamMethod);
    end
    
    varargout{1} = Def;
    varargout{2} = Mom;
         
end


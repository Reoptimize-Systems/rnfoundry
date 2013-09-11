function [g, actForce] = airgapclosure_PMSM(design, E, IVars, IMethod, sections, mode, forceRatio, initDef, alphab)
% airgapclosure_PMSM: a function for calculaing the variation in the
% air-gap due to the deflection of the supporting structure undergoing
% stresses due to the magnetic fields in the machine
%
% Input:
%
%   design - Standard PMSM design structure
%
%   E - (1 x 2) Vector containing Young's modulus of the beam materials.
%       The first value is Young Modulus of the I-beam support on the
%       stator(s) the second is for the translator.
%
%   IVars - (1 x 2) cell array of values, each of which contains a matrix
%           of values for calculating the second moment of area of the
%           stator and translator sections repectively, the the values in
%           the second cell are used for both stator sections if a double
%           sided machine is being investigated.
%
%           The contents of the first cell of IVars is determined by the
%           type of beam to be used to support the structure, as specified
%           in the IMethod field (below). 
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
%           IVars{2}(1,2): Wp, the pole width
%           IVars{2}(1,3): hs, the tooth height
%           IVars{2}(1,4): bt, the tooth width
%
%   IMethod - Defines the type of beam section used when calculating the
%             area moment of inertia. See 'MomentOfInertiaY1' for details.
%
%   sections - the number of sections into which the beam supporting the
%              magnets on the field is to be split into for the purposes of
%              calculating the variation in force with deflection
%
%   mode - a scalar value, if 0, values are calculated for a single-sided
%          machine, if 1, a double sided machine is analysed.
%
%   initDef - either (n x 1) vector or (n x sections) matrix. This contains
%             values of initial deflections in the machine parts. I.e. if
%             single sided, this is either a (2 x 1) vector or (2 x
%             sections) matrix. If vector, each value is a uniform
%             deflection along the span of the corresponding part. If
%             matrix, this is the initial deflection in each section as
%             defined by 'sections'. In either machine the deflections are
%             positive upwards and negative downwards. Conceptually the
%             machine is thought of lying on one side. 
%
%             If we consider a single-sided machne, the armature would lie
%             below the field such that if the two are being squeezed
%             toward one-another the deflection in the field is negative,
%             and the deflection in the armature is positive. The initial
%             deflections in the field are in the first row of initDef, and
%             the initial deflections in the armature in the second row.
%
%             If we consider a double-sided machne, the initial delfections
%             in the top field are in row 1, the deflections in the
%             armature are in row 2, and the deflections in the bottom
%             armature are in row 3.
%             
%   alphab - optional scalar factor for increasing the length of the beam
%            in order to incorporate the bearings. Varies from 1 upwards
%            where 1 is a beam fixed at the point where the bearings begin.
%            Is set to one if omitted
%
% Output:
%
%   g - (n x sections) matrix of new values for the air gap in the machine.
%       There will be one value of g for each of the sections being
%       considered. if a single-sided machine is being examined this will
%       be a vector, if a double sided machine, this will be a matrix with
%       two rows, one for each side of the machine.

    if nargin < 6
        mode = 1;
        initDef = [0; 0; 0];
        alphab = 1;
    else
        if mode ~= 1 && mode ~= 0
           error('mode must be 0 or 1 for single-sided or double-sided respectively') 
        end
        
        if nargin < 7
            initDef = [-0.2/1000; 0.4/1000; 0.2/1000];
            alphab = 1;
        elseif nargin < 8
            alphab = 1;
        end
    end
    
    if isempty(initDef)
        initDef = [-0.2/1000; 0.4/1000; 0.2/1000];
    end
    
    % Beam length will be stack length, plus some extra room added later
    % for bearings
    totalLength = design.ls;

    % initial Air Gap
    g = design.g;
    
    defThresh = g * 0.005;
    
    % break translator into sections
    if alphab > 1
        % add some extra space at either end of beam
        aL = [0, ((alphab-1)/2)*totalLength:totalLength/sections:((alphab-1)/2)*totalLength + totalLength, (alphab-1)*totalLength + totalLength];
        % we will want deflections in the middle of the sections so set a value z to
        % these postions
        z = [(aL(2)-aL(1))/2, (totalLength/(2*sections)) + aL(2:end-2), aL(end-1) + (aL(2)-aL(1))/2];
    else
        % beam is just the length of the stack length
        aL = 0:totalLength/sections:totalLength;
        % we will want deflections in the middle of the sections so set a value z to
        % these postions
        z = (totalLength/(2*sections)) + aL(1:end-1);
    end

    % Make sure we run once
    firstRun = true;
    
    if mode == 0
        
        % Initialise the deflections
        %Def = zeros(2,size(z,2));
        Def = ones(1,size(z,2)) .* initDef(1,:);
        Def(2,:) = ones(1,size(z,2)) .* initDef(2,:);
        
        if alphab > 1
            % get g for sections
            g = repmat(g,mode+1,size(z,2)-2) + Def(1,2:end-1) - Def(2,2:end-1);
            % Initialise section force to 0
            sectionForces = zeros(2,size(z,2)-2);  
        else
            % get g for sections
            g = repmat(g,mode+1,size(z,2)) + Def(1,:) - Def(2,:);
            % Initialise section force to 0
            sectionForces = zeros(2,size(z,2));
        end
        
        actForce = zeros(size(sectionForces));
        
        while (min(min(g)) > 0 && max(max(abs(Def))) > defThresh) || firstRun

            firstRun = false;
            
            % get new forces on sections for single sided machine
            for i = 1:size(g,2)

                % The following assumes a polynomial has been fitted to
                % the force per unit area of translator surface versus
                % the air gap
                Force = polyvaln(design.p_gforce, g(1,i));
                % Get the total force over the surface of a pole
                Force = Force * design.Wp * design.ls;
                % From this calculate the total force in each section
                Force = Force / sections;
                % Store the actual forces magnitudes
                actForce(1,i) = Force;
                % Store the actual forces magnitudes
                actForce(2,i) = -Force;
                
                % Get the increase in force due to the deflections

                sectionForces(1,i) = Force - sectionForces(1,i);

                sectionForces(2,i) = -Force - sectionForces(2,i);

            end
            
            % New air gap is reduced by deflection on both sides (opposite
            % magnets)
            if alphab > 0            
                % get new incremental increase in deflection of sections due to the
                % increase in force
                Def = beamdef_PMSM(IVars, aL(end), aL, z, [[0;0], sectionForces, [0;0]], E, IMethod, forceRatio);
                g = g - sum(abs(Def(2:end-1)));
            else
                Def = beamdef_PMSM(IVars, aL(end), aL, z, sectionForces, E, IMethod, forceRatio);
                g = g - sum(abs(Def));
            end
            
        end
    
    else
        
        % Initialise the deflections
        Def = zeros(3,size(z,2));
        Def(1,:) = ones(1,size(z,2)) .* initDef(1,:);
        Def(2,:) = ones(1,size(z,2)) .* initDef(2,:);
        Def(3,:) = ones(1,size(z,2)) .* initDef(3,:);
        
        if alphab > 1
            g = repmat(g,mode+1,size(z,2)-2);
            g(1,:) = g(1,:) + Def(1,2:end-1) - Def(2,2:end-1);
            g(2,:) = g(2,:) + Def(2,2:end-1) - Def(3,2:end-1); 
            sectionForces = zeros(3,size(z,2)-2);
        else
            g = repmat(g,mode+1,size(z,2));
            g(1,:) = g(1,:) + Def(1,:) - Def(2,:);
            g(2,:) = g(2,:) + Def(2,:) - Def(3,:);
            sectionForces = zeros(3,size(z,2));
        end
        
        actForce = zeros(size(sectionForces));

        while (min(min(g)) > 0 && max(max(abs(Def))) > defThresh) || firstRun == 1

            firstRun = false;
            
            % get new forces on sections for double sided machine  
            for j = 1:size(g,1)

                for i = 1:size(g,2)
                    % The following assumes a polynomial has been fitted to
                    % the force per unit area of translator surface versus
                    % the air gap
                    Force = polyvaln(design.p_gforce, g(j,i));
                    % Get the total force over the surface of a pole
                    Force = Force * design.Wp * design.ls;
                    % From this calculate the total force in each section
                    Force = Force / sections;
                    
                    % Get the increase in force due to the deflections
                    if j == 1 
                                                       
                        % Store the actual forces magnitudes
                        actForce(j,i) = Force;
                       
                        sectionForces(j,i) = Force - sectionForces(j,i);
                        
                    elseif j == 2
                                        
                        % Store the actual forces magnitudes
                        actForce(3,i) = Force;
                        
                        sectionForces(3,i) = -Force - sectionForces(3,i);
                        
                        % Store the net forces on armature
                        actForce(2,i) = actForce(1,i) - actForce(3,i);
                        
                        sectionForces(2,i) = sectionForces(1,i) + sectionForces(3,i);
                        
                    end

                end

            end  
            
            if alphab > 1
                % get new incremental increase in deflection of sections due to the
                % increase in force
                Def = beamdef_PMSM(IVars, aL(end), aL, z, [[0;0;0] sectionForces [0;0;0]], E, IMethod, forceRatio);
                % New air gap is reduced by deflection on both sides
                g = g + [Def(2,2:end-1) - Def(1,2:end-1); Def(3,2:end-1) - Def(2,2:end-1)];                
            else
                % get new incremental increase in deflection of sections due to the
                % increase in force
                Def = beamdef_PMSM(IVars, totalLength, aL, z, sectionForces, E, IMethod, forceRatio);
                % New air gap is reduced by deflection on both sides
                g = g + [Def(2,:) - Def(1,:); Def(3,:) - Def(2,:)];
            end
            
        end
        
    end
    % new air gaps
end

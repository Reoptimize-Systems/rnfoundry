function [g, actForce] = airgapclosure_ACPMSM(design, E, IVars, sections, IMethod, forceRatio, initDef, alphab)
% function: airgapclosure_ACPMSM
%
% Input:
%
%   E - young's modulus of the field support material
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
% design: astructure containing the following members describing the
% machine:
%
%   dgVlm - distance from back of magnet to machine centre versus the
%           magnet thickness ratio
%
%   bpVlm - magnet width to magnet depth ratio
%
%   taupVbp - pole pitch to magnet width ratio
%
%   lsVbp - machine depth to magnet width ratio
%
%   dbiVlm - field back iron thickness to magnet thickness ratio
%
%   lm - The magnet width in the same direction as the force (i.e. it's
%        thickness)
%
%   sections - the number of sections into which the beam supporting the
%              magnets on the field is to be split into for the purposes of
%              calculating the variation in force with deflection
%
% Output:
%
%   dg - (2 x sections) vector of new values for the air gaps in the machine
%        (half the distance between the two sided of the field as measured
%        from the back of the magnets in the case of the ACPMSM). There
%        will be one value of dg for each of the sections being considered.
%  
    if nargin < 9
        IMethod = '1.6';
    end

    if nargin < 8
        defThresh = design.g * 0.005;
    end
    
    if nargin < 7 || isempty(initDef)
        initDef = [0; 0];
    end

    % First calculate the basic dimensions from the ratios

    % Tooth width
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

    % Initialise the deflections
    %Def = zeros(2,size(z,2));
    Def = ones(1,size(z,2)) .* initDef(1,:);
    Def(2,:) = ones(1,size(z,2)) .* initDef(2,:);

    if alphab > 1
        % get g for sections
        g = repmat(g,1,size(z,2)-2) + Def(1,2:end-1) - Def(2,2:end-1);
        % Initialise section force to 0
        sectionForces = zeros(2,size(z,2)-2);
    else
        % get g for sections
        g = repmat(g,1,size(z,2)) + Def(1,:) - Def(2,:);
        % Initialise section force to 0
        sectionForces = zeros(2,size(z,2));
    end

    actForce = zeros(size(g));
    
    ACPMSMdir = fileparts(which('ratios2dimensions_ACPMSM'));
    load(fullfile(ACPMSMdir, 'Polynomials', 'ClosingForcePoly_ACPMSM.mat'));
    
    % Make sure we run once
    firstRun = true;
    
    while (min(min(g)) > 0 && max(max(abs(Def))) > defThresh) || firstRun

        firstRun = false;

        % get new forces on sections for single sided machine
        for i = 1:size(g,2)
            for j = 1:size(Def,1)
                
                % recalculate machine parameters
                dgVlm = ((design.Hc/2) + g(1,i)) / design.lm;

                bpVlm = design.bp / design.lm;
                taupVbp = 1 / design.bpVTaup;
                lsVbp = design.ls / design.bp;
                
                % Get the normal forces in the machine
                Force = gapclosingforce_ACPMSM(dgVlm, bpVlm, taupVbp, ...
                                         design.lsVbp, design.dbiVlm, ...
                                         design.lm, pClosing);
                
%                 % Adjust for beam of different width to pole width
%                 Force = Force .* forceRatio;

                % Get the total force over the surface of a pole
                % Force = Force * design.Taup * design.ls;
                
                % From this calculate the total force in each section
                Force = Force / sections;
                
                % Store the actual forces magnitudes
                actForce(j,i) = Force;
                
                % Get the increase in force due to the deflections
                if j == 1
                    sectionForces(j,i) = Force - sectionForces(j,i);
                elseif j == 2
                    sectionForces(j,i) = -Force - sectionForces(j,i);
                end
            end
        end

        % New air gap is reduced by deflection on both sides (opposite
        % magnets)
        if alphab > 1
            % get new incremental increase in deflection of sections due to the
            % increase in force
            Def = beamdef_ACPMSM(IVars, aL(end), aL, z, [[0;0], sectionForces, [0;0]], E, IMethod, forceRatio);
            g = g - sum(abs(Def(2:end-1)));
        else
            Def = beamdef_ACPMSM(IVars, aL(end), aL, z, sectionForces, E, IMethod, forceRatio);
            g = g - sum(abs(Def));
        end

    end

end

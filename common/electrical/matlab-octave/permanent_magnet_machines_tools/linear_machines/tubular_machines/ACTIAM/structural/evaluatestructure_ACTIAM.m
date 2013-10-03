function [PorF, z, actForce, maxDef, netResultantForce, P, design] = ...
                evaluatestructure_ACTIAM(design, options, E, netResultantForce, P)
% evaluatestructure_ACTIAM: evaluates the structure of the slotless tubular
% permanent magnet machine and reports the stress in the coil and
% deflections in the field and armature
%
% Arguments: (input)
%
%
%   E - (1 x 3) Vector containing Young's modulus of the beam materials.
%       The first value is Young Modulus of the shaft, the second the the
%       young's modulus of the coil and the third the young's modulus of
%       the outer sheath [shaft coil sheath]
%
%   Sr - the number of sections the coil is split into for evaluation
%        radially
%
%   Sz - the number of sections the coil is split into for evaluation
%        axially
%
%   avAxBVals - optional vector of the average value of the flux density in
%               the coil sections described by the Sr and Sz values at the
%               postition of maximum stress on the coil (Wp / 8 if the
%               voltage and current are in phase)
%
% Output:
%
%   PorF - (1 x 4) vector containing information about the design
%          assessment. The first value is the maximum stress in the coil,
%          the second is the maximum deflection occuring in the field, the
%          third the max deflection in the armature and the fourth, the
%          maximum total closure in the air gap. Note that this can be
%          greater than the air-gap if the first pass of the deflection of
%          the deflection algorithm results in a total deflection larger
%          than the air-gap
%
%   z - the locations along the field and armature at which the forces and
%       delections were evaluated
%
%   actForce - the total force acting between the points in z at the point
%              at which the algorithm terminated
%
%   maxDef - the criteria used to evaluate when the maximum allowed
%            deflection was breached
%
%   avAxBVals - The average value of the flux density in the coil sections
%               described by the Sr and Sz values, this will have been
%               calculated in not provided in the input, but otherwise is
%               identical to the values passed in
%    

    % perform common setup tasks for tubular machine structural evaluation
%     [design, PorF, z, peakI, peakJz, Force, avFro, avFri, Def, fFVars, ...
%         aFVars, suppWperm, polesPerSection, a1, a2, fstart, fend] = evaluatestructure_TM(design, options);
    [design, PorF, z, peakI, peakJz, Def, fFVars, aFVars, suppWperm, polesPerSection, a1, a2, fstart, fend] = evaluatestructure_TM(design, options);
    % obtain the reaction force acting at full speed for various
    % displacement values
    xshift = linspace(0, design.g - design.FEMMTol, 5);
    
    if nargin < 4
        for i = 1:size(xshift,2)
            [ReactionForce, netResultantForce(i)] = resultantforce_ACTIAM(design, design.Wp/2, peakJz, xshift(i));
            if i == 1
                P = ReactionForce;
            end
        end
        P = P * polesPerSection * options.StructSections;
    end
    
    %           IVars(1,1): Rso, outer radius of shaft support
    %           IVars(1,2): Rsi, inner radius of shaft support
    %
    %           IVars(2,1): Ro, outer coil radius
    %           IVars(2,2): Ri, inner coil radius
    %
    %           IVars(3,1): Ra, outer sheath radius
    %           IVars(3,2): Ro, inner sheath radius
    IVars = [design.Rso, design.Rsi; ...
             design.Ro, (design.Rm+design.g); ...
             design.Ra, design.Ro];
    
    % Recalculate the initial deflection taking the weight into account.
    % Here we take the absolute deflection as we want the worst case
    % scenario for the initial deflection due to weight and tolerances
    Def = Def + abs(beamdef_TM(IVars, design.supportLengths, ...
        design.totalLength, z(z>=a1 & z<=a2)-a1, P, fFVars, ...
        aFVars, z, E, suppWperm));
    
    % Calculate x, the total of the airgap closure for each z position,
    % this is the deflection from centre of the field at any position. 
    x = sum(Def,1);
    
    % We will test until g fully closes in order to compare how badly
    % different configurations fare, i.e so we can determine by how much
    % the desired maximum deflection was exceeded up to full air-gap
    % closure
    maxDef = 0.999 * design.g;
    
    % Initialise a location to store the actual total net force between the
    % field and armature in each section
    actForce = zeros(1, size(x,2)-1);
    
    DefBreakDown = Def;
    
    defThresh = design.g/1000;
    
    firstrun = true;
    
    while firstrun || (max(abs(x)) < maxDef && max(max(abs(Def))) < maxDef && max(max(abs(Def))) > defThresh)
       
        firstrun = false;

        % get new forces on sections for tubular machine  
        i = 1;
        
        for j = fstart:fend

            % Determine the total forces in the section from the per-pole
            % force
            SecNetResultantForce = interp1(xshift, netResultantForce, x(j), 'spline') * polesPerSection;

            % Store the actual forces magnitudes. This is the total force
            % acting between each section
            actForce(j) = actForce(j) + SecNetResultantForce;

            % Get the incremental increase in force due to the deflections
            fFVars(i) = SecNetResultantForce;
            % The forces on the armature will be identical to those on the
            % field but in the opposite direction
            aFVars(i) = -fFVars(i);

            i = i + 1;

        end

        % New air gap is reduced by deflection on both sides         
        % get new incremental increase in deflection of sections due to the
        % increase in force. We do not include the support lengths as we
        % are considering only the incremental forces
        Def = abs(beamdef_TM(IVars, design.supportLengths, design.totalLength, z(z>=a1 & z<=a2)-a1, P, fFVars, aFVars, z, E));
        % Store the deflection so we can apportion the point of failure,
        % i.e. field or armature
        DefBreakDown = DefBreakDown + Def;
        % Subtract the sum of the deflections in each part to the deflection
        % from centre at each position in x
        x = x + sum(Def,1);

    end

    % Pass out the difference between the maximum deflections in each
    % part and the allowed deflection, a higher function can determine
    % the appropriate course of action
    [C,I] =  max(abs(DefBreakDown(1,:)));
    PorF(2) = DefBreakDown(1,I);
    [C,I] =  max(abs(DefBreakDown(2,:)));
    PorF(3) = DefBreakDown(2,I);
    [C,I] =  max(abs(x));
    PorF(4) = abs(x(I));
    
end
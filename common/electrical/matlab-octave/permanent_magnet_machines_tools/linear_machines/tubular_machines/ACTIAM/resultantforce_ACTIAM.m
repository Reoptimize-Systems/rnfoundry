function varargout = resultantforce_ACTIAM(varargin)
% ResultantForce: This function calculates the axial reaction force from 
% the translator, the net radial force on the translator for a given
% deflection, i.e., the air-gap closing force, the total radial stress on 
% the coils and other quantities of interest
%
% Syntax:
%         
% [ReactionForce] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...,ResultantForce] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...,TotalStress] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...,TotalCoilStress] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...,TotalSheathStress] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...,ResultantCoilForce] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...,ResultantSheathForce] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x)
% 
% [...] = resultantforce_ACTIAM(Rm, WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, N, g, pos, I, x, Polynomials)
% 
% [...] = resultantforce_ACTIAM(design, pos, Jz, x)
% 
% Arguments: (input)
%
%   Rm - Translator radius in m
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Rm/Ro Ratio for machine to be evaluated, in
%   order to define the coil height
%
%   cwVWp - scalar value of Ratio of coil width to pole width (should 
%   normally be 1/3 but in the interests of generalisation will make 
%   other values up to ch/Wp = 1 possible) for machine to be evaluated
%
%   N - total number of turns in coil
%
%   g - air-gap distance in m
%
%   Rm - Radius of translator in m
%
%   pos - position of centre of coil relative to centre of a
%   steel piece in m
%
%   I - current developed at position given in pos
%
%   x - deflection of the translator from centre in m
%
% Arguments: (output)
%
%   ReactionForce - The axial reaction force from the translator in
%                   Newtons
%
%   ResultantStress - The net air-gap closing force between the translator
%                     and the stator, will be zero with no deflection from
%                     centre.
%
%   TotalStress - This is the total stress force acting on the coil and
%                 sheath
%
%   TotalCoilStress - The total stress in the coil only
%
%   TotalSheathStress - The total stress acting on the sheath only
%
%   ResultantCoilStress - The net force on the coil only in the x direction
%
%   ResultantSheathStress - The net force on the sheath only in the x
%                           direction
%
% Copyright 2009 Richard Crozier and The Institute For Energy Systems at
% The University of Edinburgh

    k = 10;
    deg = 360/k;
    delA = deg * (pi/180);

    % Initialise counter to one
    n = 1;

    load('SheathForcePoly_IA.mat');
    
    if nargin >= 12

        % Preallocate force matrix for speed
        Force = zeros(k,6);
    
        if nargin < 13
            load('BPolynomials_IA.mat');
        else
            Polynomials = varargin{13};
        end
        
        Rm = varargin{1};
        WmVWp = varargin{2};
        WpVRm = varargin{3};
        RoVRm = varargin{4}; 
        RaVRo = varargin{5};
        RsoVRm = varargin{6}; 
        WcVWp = varargin{7}; 
        N = varargin{8}; 
        g = varargin{9}; 
        pos = varargin{10}; 
        I = varargin{11}; 
        x = varargin{12}; 
        
        Sr = 8;
        Sz = 8;
        Wp = Rm * WpVRm;
        Ri = g + Rm;

        % Calculate the winding height
        ch = (RoVRm * Rm) - Rm - g;
    
        for A = 0:delA:(2*pi - delA)

            % R2 is the change in the size of the airgap at any given position
            % around the circumference when deflected and is calculated using
            % simple trigonometry
            R2(n,1) = ((x + Rm*cos(A))^2 + (Rm*sin(A))^2)^0.5;

            % Determine the new value of the airgap at this position
            g(n,1) = Ri - R2(n,1);

            % Calculate a new value for RoVRm at this position around the
            % circumference
            RoVRm(n,1) = (g(n,1) + ch + Rm) / Rm;

            % Determine the average flux density in each coil section in the
            % axial and radial directions
            %         avAxBvals = ones(9*9,1);
            %         avRadBvals = avAxBvals;meanBvals_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, Sr, Sz, Rm, pos, Wp, g)
            [avAxBvals, avRadBvals] = meanBvals_ACTIAM(WmVWp, WpVRm, RoVRm(n,1), RaVRo, RsoVRm, WcVWp, Sr, Sz, Rm, pos, Wp, g(n,1), Polynomials);

            % Determine the force acting on 1/k of a coil turn for a machine
            % with this airgapGetCoilClosingForce_TM(I, avAxBvals, RoVRm, N,
            % Rm, g, div, avFinfo)
            Force(n,1) = coilclosingforce_TM(I, avAxBvals, RoVRm(n,1), N, Rm, g(n,1), k, [Sr, Sz, Wp]);

            % Add the force due to the attraction of the ferromagnetic sheath
            Force(n,2) = sheathforcefrompolys_ACTIAM(WmVWp, WpVRm, RoVRm(n,1), RaVRo, RsoVRm, Rm, k, Polynomial);

            % Resolve the coil forces with respect to the x direction
            Force(n,3) = Force(n,1) * cos(A);

            % Resolve the coil forces with respect to the x direction
            Force(n,4) = Force(n,2) * cos(A);

            % Resolve the coil forces with respect to the x direction
            Force(n,5) = Force(n,3) + Force(n,4);

            % Determine the shear force on the coil
            Force(n,6) = coilreactionforce_TM(RoVRm(n,1), Rm, g(n,1), N, avRadBvals, I, k);

            n = n + 1;

        end

        % ReactionForce
        varargout{1} = sum(Force(:,6));
        % ResultantForce
        varargout{2} = sum(Force(:,5));
        % Total Stress
        varargout{3} = sum(Force(:,1)+Force(:,2));
        % Total Coil stress
        varargout{4} = sum(Force(:,1));
        % Total Sheath stress
        varargout{5} = sum(Force(:,5));
        % Resultant (net) coil Force
        varargout{6} = sum(Force(:,3));
        % Resultant {net) Sheath Force
        varargout{7} = sum(Force(:,4));

    elseif nargin == 4

        design = varargin{1};
        pos = varargin{2};
        Jz = varargin{3};
        x = varargin{4};

        % Initialise counter to one
        n = 1;
        
        % Preallocate force matrix for speed
        Force = zeros(k,5);

        pos = pos ./ design.Wp;

        for A = 0:delA:(2*pi - delA)

            % R2 is the change in the size of the airgap at any given position
            % around the circumference when deflected and is calculated using
            % simple trigonometry
            R2 = ((x + design.Rm*cos(A))^2 + (design.Rm*sin(A))^2)^0.5;

            % Determine the force acting on 1/k of a coil turn for a machine
            % with this airgap
            Force(n,1) = coilclosingforce_TM(design, Jz, R2 - design.Rm, k);

            g = design.Ri - R2;

            % Calculate a new value for RoVRm at this position around the
            % circumference
            RoVRm = (g + design.Hc + design.Rm) / design.Rm;
            
            % Add the force due to the attraction of the ferromagnetic
            % sheath 
            Force(n,3) = sheathforcefrompolys_ACTIAM(design.WmVWp, design.WpVRm, RoVRm, design.RaVRo, design.RsoVRm, design.Rm, k, Polynomial);
            
            % Resolve the coil forces with respect to the x direction
            Force(n,4) = Force(n,1) * cos(A);
            
            % Resolve the sheath forces with respect to the x direction
            Force(n,5) = Force(n,2) * cos(A);

            n = n + 1;

        end

        pos = pos + coilpos(design.Phases);

        % ReactionForce
        varargout{1} = sum(ylorentzforce(Jz, design.slm_intBx, pos, design.MTL));
        % ResultantForce
        varargout{2} = sum(Force(:,4) + Force(:,5));
        % TotalStress
        varargout{3} = sum(Force(:,1) + Force(:,2));
        % Total Coil stress
        varargout{4} = sum(Force(:,1));
        % Total Sheath stress
        varargout{5} = sum(Force(:,2));
        % Resultant (net) coil Force
        varargout{6} = sum(Force(:,4));
        % Resultant {net) Sheath Force
        varargout{7} = sum(Force(:,5));
        
    else

        error('Incorrect number of arguments')
    end

end


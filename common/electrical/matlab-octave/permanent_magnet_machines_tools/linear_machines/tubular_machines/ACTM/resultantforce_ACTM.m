function varargout = resultantforce_ACTM(varargin)
% ResultantForce: This function calculates the axial reaction force from 
% the translator, the net radial force on the translator for a given
% deflection, i.e., the air-gap closing force, the total radial stress on 
% the coils 
%
% Syntax
% 
% [ReactionForce] = resultantforce_ACTM(Rm, WmVWp, WpVRm, RoVRm, RsoVRm, WcVWp, N, g, pos, I, x, Polynomials)
% 
% [..., ResultantForce] = resultantforce_ACTM(Rm, WmVWp, WpVRm, RoVRm, RsoVRm, WcVWp, N, g, pos, I, x, Polynomials)
% 
% [..., TotalStress] = resultantforce_ACTM(Rm, WmVWp, WpVRm, RoVRm, RsoVRm, WcVWp, N, g, pos, I, x, Polynomials)
% 
% [ReactionForce] = resultantforce_ACTM(design, pos, Jz, x)
% 
% [..., ResultantForce] = resultantforce_ACTM(design, pos, Jz, x)
% 
% [..., TotalStress] = resultantforce_ACTM(design, pos, Jz, x)
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
%           order to define the coil height
%
%   WcVWp - scalar value of Ratio of coil width to pole width (should 
%           normally be 1/3 but in the interests of generalisation will
%           make other values up to ch/Wp = 1 possible) for machine to be
%           evaluated
%
%   N - total number of turns in coil
%
%   g - air-gap distance in m
%
%   Rm - Radius of translator in m
%
%   pos - position of centre of coil relative to centre of a
%         steel piece in m
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
%   TotalStress - This is the total stress force acting on the coil
%
% Copyright 2007 Richard Crozier and The Institute For Energy Systems at
% The University of Edinburgh
        
    k = 10;
    deg = 360/k;
    delA = deg * (pi/180);
        
    if nargin >= 11
        
        Rm = varargin{1}; 
        WmVWp = varargin{2}; 
        WpVRm = varargin{3}; 
        RoVRm = varargin{4}; 
        RsoVRm = varargin{5}; 
        WcVWp = varargin{6}; 
        N = varargin{7}; 
        g = varargin{8}; 
        pos = varargin{9}; 
        I = varargin{10}; 
        x = varargin{11}; 
        
        Sr = 8;
        Sz = 8;
        Wp = Rm * WpVRm;
        Ri = g + Rm;

        % Calculate the winding height
        ch = (RoVRm * Rm) - Rm - g;

        % Initialise counter to one
        n = 1;

        % Preallocate force matrix for speed
        Force = zeros(k,3);

        if nargin < 12
            load('BPolynomials_AC.mat');
        else
            Polynomials = varargin{12};
        end

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
            %         avRadBvals = avAxBvals;
            [avAxBvals, avRadBvals] = GetAvBvals_ACTM(WmVWp, WpVRm, RoVRm, RsoVRm, WcVWp, Sr, Sz, Rm, pos, Wp, g, Polynomials);

            % Determine the force acting on 1/k of a coil turn for a machine
            % with this airgap
            Force(n,1) = coilgapclosingforce_TM(I, avAxBvals, RoVRm(n,1), N, Rm, g(n,1), k);

            % Resolve the force with respect to the x direction
            Force(n,2) = Force(n,1) * cos(A);

            % Determine the shear force on the coil
            Force(n,3) = coilreactionforce_TM(RoVRm(n,1), Rm, g(n,1), N, avRadBvals, I, k);

            n = n + 1;

        end

        % ReactionForce
        varargout{1} = sum(Force(:,3));
        % ResultantForce
        varargout{2} = sum(Force(:,2));
        % TotalStress
        varargout{3} = sum(Force(:,1));

    elseif nargin == 4
        
        design = varargin{1};
        pos = varargin{2};
        Jz = varargin{3};
        x = varargin{4};

        % Preallocate force matrix for speed
        Force = zeros(k,2);
        % get the angles of the coil section around the circumference
        A = (0:delA:(2*pi - delA))';
        % R2 is the change in the size of the airgap at any given position
        % around the circumference when deflected and is calculated using
        % simple trigonometry
        R2 = sqrt((x + design.Rm.*cos(A)).^2 + (design.Rm.*sin(A)).^2);
        % Determine the force acting on 1/k of a coil turn for a machine
        % with this airgap
        Force(:,1) = coilclosingforce_TM(design, Jz, R2 - design.Rm, k);
        % Resolve the force with respect to the x direction
        Force(:,2) = Force(:,1) .* cos(A);
        
        pos = pos ./ design.Wp;

        pos = pos + coilpos(design.Phases);
        
        % ReactionForce
        varargout{1} = sum(ylorentzforce(Jz, design.slm_intBx, pos, design.MTL));
        % ResultantForce 
        varargout{2} = sum(Force(:,2));
        % TotalStress
        varargout{3} = sum(Force(:,1));
        
    else

        error('Incorrect number of arguments')
    end

end


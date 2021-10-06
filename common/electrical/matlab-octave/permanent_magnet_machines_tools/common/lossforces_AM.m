function [Ftot, Firon, Feddy] = lossforces_AM(design, simoptions, xR, vR)
% calculate the forces due to losses in the components of an electrical
% machine
%
% Syntax
%
% Fadd = lossforces_AM(design, simoptions, xR, vR)
%
% Description
%
% lossforces_AM calculates the force due to losses in an electrical machine
% due to iron losses and eddy current losses in the windings. It requires
% that certain precalculated information is supplied in the form of several
% slm objects. The inputs design and simoptions are structures containing
% information necessary for the calculation of these losses. The simoptions
% structure should contain the field:
%
%   NoOfMachines - a scalar value of the the number of machines connected
%     together.
%
% The design structure must contain the following fields:
%
%   PowerPoles - a scalar value of the number of Poles involved in the iron
%     loss in the machine
%
%   NStages - the number of stages in a multistage machine
%
%   NCoilsPerPhase - the number of coils per phase in the machine
%
%   slm_eddysfdpart - this field must contain an slm object fitted to the 
%     part calculation of winding eddy current losses in a coil due to an
%     externally applied time-varying field. The resulting calculation of
%     losses is based on the Squared Field Derivative (SFD) method
%     presented in [1]. The SFD method uses the formula:
%                              
%         pi * l * N * d_c^4   
%     P = ------------------ . < ( dB / dt )^2 > 
%              64 rho_c         
%                              
%     Where l is the length of a strand in a winding in the changing field,
%     N is the product of the number of strands and number of turns, d_c is
%     the diameter of each strand, rho_c is the resistivity of the wire
%     material, and < ( dB / dt )^2 >  is the spacial average of the value
%     of ( dB / dt )^2 over the winding cross-section, to calculate the
%     instantaneous power losses in the coil. lossforces_Am exploits the
%     fact that for electrical machines, the value of dB / dt at any time
%     is closely related to the motion of the machine such that it is
%     equivalent to (dB / dx)*(dx/dt) = (dB / dx)*(v). Therefore, for a
%     machine winding the losses can be calculated using the formula:
%
%               pi * l * N * d_c^4   
%     P = v^2 * ------------------ . < ( dB / dx )^2 > 
%                   64 rho_c         
%                        
%     The field slm_eddysfdpart is expected to contain a periodic SLM
%     object fitted to the value of the part of this expression to the left
%     of v^2 at positions normalised to the pole width (i.e. fitted to
%     values of the expression between normalised x values between 0 and
%     2), such that multiplying the value returned by evaluating the slm by
%     v^2 yields the instantaneous losses in a winding coil due to eddy
%     currents in the conductors.
%
%   CoreLossSLMs - this is a structure containing three fields which are
%     each SLM objects fitted to functions used to determine the per-pole
%     iron losses in the machine. These losses are calculated using the
%     method described in [2]
%
%
% The direction of the calculated forces due to the losses is negative when
% the velocity is positive and vice versa.
%
% [1] C. R. Sullivan, "Computationally efficient winding loss calculation
% with multiple windings, arbitrary waveforms, and two-dimensional or
% three-dimensional field geometry," IEEE Transactions on Power
% Electronics, vol. 16, no. 1, pp. 142--150, 2001.
%
% [2] D. Lin, P. Zhou, W. N. Fu, Z. Badics, and Z. J. Cendes, "A Dynamic
% Core Loss Model for Soft Ferromagnetic and Power Ferrite Materials in
% Transient Finite Element Analysis," IEEE Transactions on Magnetics, vol.
% 40, no. 2, pp. 1318--1321, Mar. 2004.
%
%
% See also: softferrolossrectregionvarxpartcalc.m
%

% Copyright Richard Crozier 2012-2013
    
    Pironloss = corelossfcn(abs(vR), xR, design.CoreLossSLMs.hxslm, design.CoreLossSLMs.cxslm, design.CoreLossSLMs.exslm) ...
                    * design.PowerPoles * design.NStages * simoptions.NoOfMachines;

    Pexteddy = sum(periodicslmeval(xR + design.CoilPositions .* design.PoleWidth, design.slm_eddysfdpart, 0, false)) ...
                    .* realpow(abs(vR),2) * design.NCoilsPerPhase * design.NStages * simoptions.NoOfMachines;

	Firon = power2force(vR, Pironloss);
    
    Feddy = power2force(vR, Pexteddy);
    
    Ftot = Firon + Feddy;
    
end


function P = corelossfcn(v, x, hxslm, cxslm, exslm)
% evaluates the soft ferrite losses in the core from partially
% precalculated loss functions
%
% Syntax
%
% P = corelossfcn(v, x, hxslm, cxslm, exslm)
%
% Input
%
%   v - relative velocity of the machine parts
%
%   x - relative position of the machine parts
%
%   hxslm - slm object fitted to the hysteresis loss function over two
%     Poles
%
%   cxslm - slm object fitted to the eddy current loss function over two
%     Poles
%
%   exslm - slm object fitted to the excess loss function over two Poles
%
%
% Output
%
%   P - the instantaneous power loss in the iron core
%
    
% Copyright Richard Crozier 2012-2013

    P = abs(v .* periodicslmeval(x, hxslm, 0, false) ...
        + realpow(v,2) .* periodicslmeval(x, cxslm, 0, false) ...
        + realpow(sqrt(v),3) .* periodicslmeval(x, exslm, 0, false));
    
    % smooth the function around zero velocity for the ode solvers. Without
    % introducing a gradient below a threshold the forces calculated from
    % the losses flip from positive to negative about zero velocity. This
    % causes problems for the ode solvers which attempt to take smaller and
    % smaller time steps to maintain the acuracy of the solution.
    vthresh = 0.01;
    % the function v^2 / vthresh^2 is a quadratic which rises from zero to
    % 1 between v = 0 and v = vthresh. 
    P(v < vthresh) = P(v < vthresh) .* v(v < vthresh).^2 ./ (vthresh^2); 
    
end
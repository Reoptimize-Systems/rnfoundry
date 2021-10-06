function [dpsidxF, EMF, Feff, ForceVec, xE, vE, unitv, design] = machineodesim_linear(design, simoptions, Icoils, xBh, vBh, xBs, vBs)
% determines the EMF, rate of change in flux linkage with relative
% armature/translator position, Shear Force and the velocities of each part
% of the machine
%
% Syntax
%
% [dpsidxF, EMF, Feff, ForceVec, xE, vE] = ...
%       machineodesim_linear(design, simoptions, Icoils, xBh, xBs, vBh, vBs)
%
% Input
%
% design and simoptions are structures which include certain fields
% necessary for the machine simulation. design must contain at least the
% following fields:
%   
%    FieldDirection: a scalar value of either 1 or -1 indicating the
%      direction of the magnetic field displacement relative to the
%      direction of the prime mover displacement. If 1 this means the
%      magnetic field moves in the same direction as the translator, i.e.
%      the magnets are mounted on the translator and the coils are
%      stationary. If -1, the opposite is true, and the magnetic field is
%      moving in the opposite direction.
%
%    PoleWidth: scalar value of the physical width of one magnetic pole of
%      the machine.
%
%    CoilPositions: vector of values contianing the relative positions of
%      the coils in a phase, normalised to the pole width. For example, 
%      design.CoilPositions = [0, 2/3, 4/3]. Some typical normalised coil
%      spacings are provided by coilpos.m for different phase numbers.
%
%    slm_psidot: An slm object (as produced by slmengine) fitted to the
%      flux linkage in a macine coil from the position of maximum flux
%      linkage to the position of minimum flux linkage against the
%      normalised displacement over one pole.
%
%    PowerPoles: scalar value of the number of active, power producing
%      Poles in the machine
%
%    TemperatureBase: scalar value of a base temperature at which a
%      conductor resistivity will be supplied, and temperature coefficient
%      allowing the temperature dependent resistance to be determined.
%
%    WireResistivityBase: scalar value containing the base value of the
%      resistivity of the conductors in the coil at the supplied base
%      temperature in TemperatureBase.
%
%    AlphaResistivity: temperature coefficient of the conductor resistivity
%      supplied in WireResistivityBase such that the actual resistivity is
%      given by:
%
%      rho = rhobase * (1 + alpha * (T - Tbase));
%
%    RDCBase: Matrix of base values of the resistance in of each phase at
%      the base temperature in TemperatureBase. This must be a [n x n]
%      matrix of values for n Phases, where the diagonal terms are the
%      resistance of each phase, and the off-diagonal terms are typically
%      zero, e.g. 
%
%      design.RDCBase = [R1, 0, 0; 0, R2, 0; 0, 0, R3]
%
%    Phases - scalar number of Phases in the machine armature.
%
%
%  simoptions - should contain at least the following fields:
%
%    tether_length: this is the starting lenght of a rigid tether from the
%      prime mover to the effector of the machine.
%
%    NoOfMachines: the number of machines linked together for the purposes
%      of force calculation.
%
%    Temperature: scalar value of the temperature of the coil conductors,
%      used to determine the temperature-dependent resistance.
%
%
% The following arguments are then:
%
% Icoils - A vector of values of the the coil current, should be the same
%   length as the number of Phases in the machine.
% 
% xBh - The displacement of the prime mover (typically a buoy) in the
%   vertical (heave) direction.
%
% xBs - The displacement of the prime mover (typically a buoy) in the
%   horizontal (surge) direction.
% 
% vBh - The velocity of the prime mover (typically a buoy) in the
%   vertical (heave) direction.
% 
% vBs - The velocity of the prime mover (typically a buoy) in the
%   horizontal (surge) direction.
%
% Output
% 
% dpsidxF -  vector of values of the derivative of the flux linkage w.r.t.
%   displacement at the current position of the machine.
% 
% EMF - vector of values of the EMF in a single coil in each phase
% 
% Feff - The force acting on the effector of the machine (the part attached
%   to the prime mover). This force is defined as positive in the direction
%   of positive displacement of the effector.
% 
% ForceVec - A vector of value of the force acting on the prime mover in
%   the surge and heave directions
% 
% xE -  The vertical displacement of the effector of the machine.
% 
% vE - The vertical velocity of the effector of the machine.
%

    % first determine the translator displacement from the tether
    % length, the change in translator vertical position will be the change
    % in distance from hawser to buoy, i.e the change in tether length
    xE = sqrt ((xBh + simoptions.BuoySim.tether_length).^2 + xBs.^2) - simoptions.BuoySim.tether_length;

    % Find unit vector in the direction pointing from hawse hole to the
    % buoy, we add a tiny length to the tether length in case it is of
    % length zero which would make calculation of the unit vector
    % impossible
    unitv = [simoptions.BuoySim.tether_length+1e-6+xBh, xBs] ...
              / norm([simoptions.BuoySim.tether_length+1e-6+xBh, xBs]);

    % Then find dot product of heave and surge velocities with unit vector
    % to get correct direction and magnitude of the velocity of the
    % magnetic field relative to the coils
    vE = dot (unitv, [vBh, vBs]);

    % Get the EMF and forces using the core machine simulation function
    [Feff, ~, EMF, dpsidxF, design] = machineodesim_AM (design, simoptions, xE, 0, vE, 0, Icoils);
    
    %  Calculating the proportion that the machine forces are producing on
    %  surge and in heave. Do this by multiplying unit vector in direction
    %  of hawser to buoy by total forces. This results in new vector:
    %
    %               FfeaVec = [Ffea_heave, Ffea_surge]
    %
    %  Which should also have forces in the right directions. The FEA
    %  forces are calculated for the coils, therefore when the translator
    %  is moving up, Ffea is also acting upwards. We must reverse their
    %  direction to get the forces acting on the translator (and buoy).
    ForceVec = Feff * unitv;

end


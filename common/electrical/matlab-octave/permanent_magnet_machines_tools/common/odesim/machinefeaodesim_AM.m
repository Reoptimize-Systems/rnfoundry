function [flux_linkage, FEF, FRE, lossinfo, design] = machinefeaodesim_AM(design, simoptions, xEF, xRE, vEF, vRE, Icoils)
% core permanent magnet machine dynamic ode simulation function
%
% Syntax
%
% [FEF, FRE, EMF, dpsidxCRTF, design] = machineodesim_AM(design, simoptions, xEF, xRE, vEF, vRE, Icoils)
%
% Description
%
% machineodesim_AM calculates the coil emfs and forces/torques on the coils
% in the phases of a permanent magnet machine.
%
% Input
%
%  design - a structure containing several fields with information about
%    the machine design. design should contain at least the following
%    fields:
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
%  simoptions - a structure containing various parameters of the
%    simulation. It should contain at least the following fields:
%
%    NoOfMachines: scalar value of the number of physically connected
%      machines in the simulation. This is used to determine the total
%      force from the linked machines.
%
%    Temperature: scalar value of the temperature of the coil conductors,
%      used to determine the temperature-dependent resistance.
%
%  xEF - position of the effector (the part of the machine directly
%    connected to the prime mover) relative to a global reference.
%
%  xRE - position of the reactor (the part of the machine NOT directly
%    connected to the prime mover) relative to a global reference.
%
%  vEF - velocity of the effector (the part of the machine directly
%    connected to the prime mover) relative to a global reference.
%
%  vRE - velocity of the reactor (the part of the machine directly
%    connected to the prime mover) relative to a global reference.
%
%  Icoils - vector of values containing the current in a coil of each phase
%    of the machine.
%
% Output
%
%  FEF - The total force applied by the effector of the machine
%
%  FRE - The total force applied to the reactor of the machine
%
%  EMF - The emf produced in a single coil in each phase of the machine
%
%  dpsidxCRTF - the derivative of the flux linkage in each phase of the
%    machine
%
%  design - the design structure, but with a field 'R' appended containing
%    the modified, temperature dependendent resistance values. 
% 
%
% See also: machineodesim_linear, machineodesim_linear_mvgarm,
%           machineodesim_linear_mvgfield
%

% Copyright Richard Crozier 2016

    % Calculate the position of the field relative to the armature. For
    % both parts positive displacement is in the same direction as positive
    % forces acting on both parts. 
    %
    posCoilRelToField = xEF - xRE;
    
    % Next we create an fea simulation of the machine with this rotor
    % position and coil currents
    [flux_linkage, FEF, lossinfo] = ...
        feval (design.FEAFluxLinkageFCN, design, simoptions, posCoilRelToField, Icoils);

    % the force on the reactor is the reverse of the effector
    FRE = -FEF;
    
    design = machineodesim_common_AM (design, simoptions, vEF - vRE);

end

            

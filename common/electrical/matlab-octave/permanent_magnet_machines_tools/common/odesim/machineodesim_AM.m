function [FEF, FRE, EMF, dpsidxCRTF, design, pos] = machineodesim_AM(design, simoptions, xEF, xRE, vEF, vRE, Icoils)
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

% Copyright Richard Crozier 2009 - 2014

    % Calculate the position of the field relative to the armature. For
    % both parts positive displacement is in the same direction as positive
    % forces acting on both parts. We normalise the position to the pole
    % width, as it is expected that the provided slm object is fitted to the
    % flux linkage of a coil with the coil displaced over one pole width
    % normalised to the domain 0 to 1.
    %
    % We also set the direction of the relative displacement or the coil to
    % the field according to the definition of the direction in the design
    % structure (design.FieldDirection) and convert to a position relative
    % to pole width. The value of design.FieldDirection should be either 1,
    % or -1. If 1 this means the magnetic field moves in the same direction
    % as the effector, e.g. the magnets are mounted on the translator and
    % the coils are stationary. Effectively this means the displacement is
    % the inverse of the effector motion. If -1, the opposite is true, and
    % the magnetic field is moving in the opposite direction. e.g. the
    % magnets are mounted on a stationary reactor and the coils on a moving
    % effector.
    xCoilRelToField = -design.FieldDirection .* (xEF - xRE) ./ design.PoleWidth;
    
    % Next we get the positions of the magnetic field relative to each coil
    % in a multi-phase block. We will use this to get the differential of
    % the flux linkage w.r.t. the magnetic field position, i.e. 
    %%
    % <latex>
    % \begin{equation}
    %   \frac{\text{d} \, \lambda}{\text{d} x_F}
    % \end{equation}
    % </latex>
    % 
%     pos = (((1:design.Phases) .* (2/design.Phases)) - 1 -
%     (1/design.Phases) + xRF)';
%     pos = (((1:design.Phases) .* (1/design.Phases + 1) )  + xCoilRelToField)';
    pos = bsxfun (@plus, design.CoilPositions(:), xCoilRelToField);
    
    % Find dpsidxR from an slm object fitted to the flux linkage versus
    % the positions 0 <= xR <= 1.0 where xR is xEF ./ design.PoleWidth, i.e.
    % the displacement of the magnetic field relative to the armature coils
    % normalised to the pole width
%     dpsidxCRTF = slmpsidot_linear(design, pos, design.PoleWidth);
    dpsidxCRTF = periodicslmeval (pos, design.slm_fluxlinkage, 1, false) ./ design.PoleWidth;
    
    % The EMF (voltage) produced in the coils is given by -ve the rate of
    % change of flux linkage with respect to the displacement of the coil
    % multiplied by the velocity of the coils magnetic field relative
    % magnetic field.
    %
    % We will first calculate the velocity of the field relative to the
    % armature
    vCoilRelToField = -design.FieldDirection .* (vEF - vRE);
    
    %%
    % <latex>
    % \begin{equation}
    %   \text{EMF} = - frac{\text{d}\,\lambda}{\text{d} x_F} * v_F
    % \end{equation}
    % </latex>
    EMF = bsxfun (@times, -dpsidxCRTF, vCoilRelToField) .* design.CoilsPerBranch;
    
    % determine the forces due to the magnets and electrical forces at
    % the relative position xR absed on the coil current and rate of change
    % of flux linkage w.r.t. xR at this point.
    %%
    % <latex>
    % \begin{equation}
    %   F_{\text{trans}} = \sum I \frac{\text{d} \, \lambda}{\text{d} x_F} 
    % \end{equation}
    % </latex>
    %
    FEF = -design.FieldDirection ...
                * sum (bsxfun (@times, Icoils, dpsidxCRTF), 1) ...
                .* design.NCoilsPerPhase ...
                .* design.NStages ...
                .* simoptions.NoOfMachines;

    % the force on the reactor is the reverse of the effector
    FRE = -FEF;
    
    % do some common sim tasks, including applying temp dependance on
    % resistivity
	design = machineodesim_common_AM (design, simoptions, vCoilRelToField);

end


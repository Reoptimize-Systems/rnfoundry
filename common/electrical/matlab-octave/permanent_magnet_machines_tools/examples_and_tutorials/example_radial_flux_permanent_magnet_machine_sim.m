%% example_radial_flux_permanent_magnet_machine_sim.m
%
% This script sets up an example simulation of a radial flux permanent
% magnet machine using the RenewNet Foundry permanent magnet machine
% toolbox and also demonstrates some other capabilities of the toolbox,
% such as the underlying common functions used by the machine simulations.
%
%
%
% A machine design is always defined using fields of a structure which we
% usually name 'design'. Throughout the toolbox this is the name used. In
% this example we will set up a simulation of a radial flux permanent
% magnet machine, and take a look at what is actually produced by the
% toolbox. 
%
% The machine we will simulate is a 3-phase, six-pole machine with an
% overlapping winding. The machine design is arbitrary and probably not a
% particularly good example of a machine design, and is purely for
% demonstrating the tools.
% 
%% Setting up the design
%
% So first we set up some design variables describing the machine in the
% structure. The first thing we do is set up the winding arrangement. The
% way this is specified is based on the Stator slot/coil/winding
% terminology presented in:
%
% J. J. Germishuizen and M. J. Kamper, "Classification of symmetrical
% non-overlapping three-phase windings," in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
%
% To summarise this the followng variables describe a winding, we don't
% need to specify all of these though, some are calculated:
%
% Poles - total number of magnetic Poles in the machine (the number of
%         magnets in permanent magnet machine)
% Phases - the number of electrical Phases in a machine
% Qs  -  total number of stator slots in all Phases combined
% Qc  -  total number of winding coils in all Phases combined
% yp - Average coil pitch as defined by (Qs/Poles)
% yd - Actual coil pitch as defined by round(yp) +/- k
% q  -  number of slots per pole and phase
% qn  -  numerator of q
% qd  -  denominator of q
% qc - number of coils per pole and phase, i.e. the ratio coils / (Poles * Phases)
% qcn  -  numerator of qc
% qcd  -  denominator of qc
%
% In addition, coils can either be single layered (only one coil side per
% slot, not overlapping in the slots), or double layered (two coils sit on
% top of each other in a slot). In each case the following relationships
% hold:
%
% Single layer
%   q = 2qc
%   Qs = 2Qc
% Double layer
%   q = qc
%   Qs = Qc
%
% With this new knowledge we can set up the winding design we will use for
% our example machine.

% Before we start, we'll get rid of anything in the workspace that might
% confuse us, i.e. any existing design and simoptions structure
clear design simoptions

% we will design a 6-pole machine, so it will have 6 magnets, and 3
% pole-pairs
design.Poles = 6;
% Choose the number of Phases, the conventional 3
design.Phases = 3;
% The desired number of layers is stored in 'CoilLayers'
design.CoilLayers = 2; 
% The type of winding is specified as a string, it can be 'overlapping' or
% 'nonoverlapping'
design.WindingType = 'overlapping';
% Next we can specify the winding by stating the ratio of coils to Poles
% and slots. This value must be a fraction object, which is created as
% below. To simplify things we will use one coil per pole and phase, but
% many other ratios are possible, see the help for "fr" (run "help fr" at
% the command prompt) for more information on using fractions objects
design.qc = fr(1,1);
% We must also specify a fill factor for the coils, this is the wire fill
% factor acheived
design.CoilFillFactor = 0.6;
% We must also specify either the number of turns, of the diameter of the
% wire used in the coils. If we specify the turns, the wire diameter is
% calculated, and vice-versa. The wire diameter is specified in the field
% "Dc" and the number of turns in the field "CoilTurns". In this case we
% will specify the number of turns
design.CoilTurns = 20;

% Now set up actual dimensions of the machine components. There are
% actually three different ways to do this. One is to specify the outer
% radius of the machine, and then the radial lengths of various components
% and other dimensions. This is the method we use here. Another way is to
% specify the radial distance from the centre of the machine of the
% components, or alternatively specify the outer radius, then all the other
% dimensions through a set of dimensionless ratios based on this first
% dimension. There are various reasons this can be more convenient.
%
% However, as stated, here we use the lengths and the outer rotor radius.
% What these dimensions are will be more obvious when we open the design
% for plotting/viewing later.

% The outer radius of the rotor back iron
design.Rbo = 0.5;
% The height of the magnets in the radial direction
design.tm = 0.005;
% The thickness of the rotor back iron
design.tbi = 0.01;
% The thickness of the stator yoke (the part which the slots sit on)
design.ty = 0.015;
% The radial height of a slot opening
design.tc = 0.03;
% The height of the slot shoe at the point it joins the slot wall
design.tsb = 0.01;
% The height of the slot shoe at the point it joins the slot wall
design.tsg = 0.003; 
% The radial air gap length between rotor and stator
design.g = 3/1000;
% magnet pitch in radians
design.thetam = (2*pi / design.Poles) * 0.8;
% coil slot opening pitch in radians ( space a coil has to fit in a slot )
design.thetac = (2*pi / design.Poles / design.Phases) * 0.85;
% the pitch in radians of gap between the shoes that partially cover the
% slot openings
design.thetasg = design.thetac * 0.8;
% stack length of the machine, the depth along the cylindrical axis of the
% machine
design.ls = 0.3;

% next we tell the simulation what type of stator/rotor arrangement we
% have. There are two options at the moment, either a stator on the inside
% "facing" outwards, i.e. with the slots pointing outwards and with a rotor
% on the outside, or we can have the opposite arrangement with an internal
% rotor and internally facing stator. 
design.StatorType = 'so';

% We must also specify some materials in the machine
design.MagnetMaterial = 'NdFeB 32 MGOe';
design.BackIronMaterial = '1117 Steel';
design.YokeMaterial = design.BackIronMaterial;
design.CoilMaterial = '36 AWG';

% for convenience, a function is provided to complete a lot of design
% parameters of the machine from the variables provided above. This
% function is completedesign_RADIAL_SLOTTED. It requires you to supply both
% a design structure and a simoptions structure. We won't worry about
% simoptions for now, this can be empty for the moment. MOre on simoptions
% later.

% create the empty simoptions structure
simoptions = struct();
% complete the design
design = completedesign_RADIAL_SLOTTED(design, simoptions);

%% Simulating the Machine
%
% The first stage in simulating the machine is to gather electromagnetic
% data using FEMM, or mfemm. The function simfun_RADIAL_SLOTTED does this
% for this type of machine. For other machine types the function name
% follows the same pattern. These functions are found in the odesim
% directory for each machine. 
%
%
% At this point we reintroduce the simoptions structure mentioned earlier.
% This structure is intended to specify various simulation options.

% GetVariableGapForce is a flag which determines whether
% simfun_RADIAL_SLOTTED does a series of force simulations at progressively
% smaller air gaps to create data on how the force varies with respect to
% this. We set this to false to save time.
simoptions.GetVariableGapForce = false;

% now run simfun_RADIAL_SLOTTED, this function performs various simulations
% of the machine
[design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions);









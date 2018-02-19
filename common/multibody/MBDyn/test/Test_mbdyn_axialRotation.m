
position = [0;0;0];
orientation = mbdyn.pre.orientmat ('eye'); 
options.AbsoluteVelocity = [0;0;0];
options.AbsoluteAngularVelocity = [0;0;0];


design.gi = 6e-3; % gap between inner rotor and fmps
design.go = 6e-3; % gap between outer rotor and fmps
design.tfmp = 64e-3; % fmp thickness
design.tmi = 4e-3; % inner rotor magnet thickness
design.tmo = 2e-3; % outer rotor magnet thickness
design.tbi = 10e-3; % inner rotor back iron thickness
design.tbo = 10e-3; % outer rotor back iron thickness
design.Rii = 100e-3;  % inner rotor, inner radius
design.Rio = design.Rii + design.tbi + design.tmi; % inner rotor, outer radius
design.Roi = design.Rio + design.gi + design.tfmp + design.go; % outer rotor, inner radius
design.Roo = design.Roi + design.tmo + design.tbo; % outer rotor, outer radius
design.ls = 2 * design.Roo; % axial length (fmp stack length)
design.Npi = 4; % inner rotor num poles
design.Npo = 16; % outer rotor num poles
design.Rs = 10e-3; % shaft radius
design.GearRatio = 10; % gear ratio (could be calculated from poles, or fea)
design.TorqueLowSpeedPeak = 1e3; % low speed peak torque


design.MassInner = 7600 * (circlearea (design.Rio) - circlearea (design.Rii)) * design.ls;
design.MassOuter = 7600 * (circlearea (design.Roo) - circlearea (design.Roi)) * design.ls;

% inertia formulas from:
%https://en.wikipedia.org/wiki/List_of_moments_of_inertia
%
% the inertia is specified in the frame of the attached node (unless
% otherwise specified), so we input the inertial matrix with the axis
% lying coincident with the x axis in the node frame.
    %
Ixx = 0.5 * design.MassInner * ( design.Rii.^2 + design.Rio.^2 );
Iyy = (1/12) .* design.MassInner .* (3 .* ( design.Rii.^2 + design.Rio.^2 ) + design.ls.^2 );
Izz = Iyy;

design.InertiaInner = diag ( [Ixx, Iyy, Izz] );

design.CoGInner = [0;0;0];

Ixx = 0.5 * design.MassOuter * ( design.Roi.^2 + design.Roo.^2 );
Iyy = (1/12) .* design.MassOuter .* (3 .* ( design.Roi.^2 + design.Roo.^2 ) + design.ls.^2 );
Izz = Iyy;

design.InertiaOuter = diag ( [Ixx, Iyy, Izz] );

design.CoGOuter = [0;0;0];

intermediate_node = mbdyn.pre.structuralNode6dof ('dynamic', ...
    'AbsolutePosition', position, ...
    'AbsoluteOrientation',  orientation, ...
    'AbsoluteVelocity', options.AbsoluteVelocity, ...
    'AbsoluteAngularVelocity', options.AbsoluteAngularVelocity);

outer_rotor_node = mbdyn.pre.structuralNode6dof ('dynamic', ...
    'AbsolutePosition', position, ...
    'AbsoluteOrientation',  orientation, ...
    'AbsoluteVelocity', options.AbsoluteVelocity, ...
    'AbsoluteAngularVelocity', options.AbsoluteAngularVelocity);

% we will make hinge always be on the x axis in the frame of the nodes
hinge_orientation = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [0;1;0], 'ib', 2, 'vecB', [0;0;1]));

nodeoffset = [0;0;0];

% formula is torque on rotor with variable being relative angle of
% revolution of inner and outer rotor

outer_rotor_drive = mbdyn.pre.multDrive ( mbdyn.pre.nodeDrive ( intermediate_node, mbdyn.pre.directDrive (), 'String', 'omega[1]' ), ...
                                          mbdyn.pre.const (design.GearRatio - 1) );


jAxRot = mbdyn.pre.axialRotation ( intermediate_node, outer_rotor_node, outer_rotor_drive, ...
                                   'Offset1', nodeoffset, ...
                                   'Offset2', nodeoffset, ...
                                   'Offset1Reference', 'node', ...
                                   'Offset2Reference', 'other node', ...
                                   'RelativeOrientation1', hinge_orientation, ...
                                   'Orientation1Reference', 'node', ...
                                   'RelativeOrientation2', hinge_orientation, ...
                                   'Orientation2Reference', 'other node' ...
                                 );


outer_rotor_body = mbdyn.pre.body ( design.MassOuter, ...
                                    design.CoGOuter, ...
                                    design.InertiaOuter, ...
                                    outer_rotor_node );

nodes = {intermediate_node, outer_rotor_node};

bodies = {outer_rotor_body};

elements = { jAxRot };



% we will make hinge always be on the x axis in the frame of the nodes
pin_orientation = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [0;1;0], 'ib', 2, 'vecB', [0;0;1]));

pinjoint = mbdyn.pre.revolutePin (nodes{1}, position, position, ...
                    'NodeRelativeOrientation', pin_orientation, ...
                    'NodeRelativeOrientationReference', 'global', ...
                    'PinOrientation', pin_orientation);
                

tstart = 0;
tend = 30;
tramp = tend/4;
tstep = 50e-3;
omega_start = 0;
omega_final = rpm2omega (10);
omega_slope = (omega_final - omega_start) ./ (tramp - tstart);

omega_ramp = mbdyn.pre.rampDrive (tstart, tramp, omega_slope, omega_start);

omega_joint = mbdyn.pre.angularVelocity (nodes{1}, [1;0;0], omega_ramp);

prb = mbdyn.pre.initialValueProblem (tstart, tend, tstep, 'ResidualTolerance', 1e-6); %, 'DerivativesTolerance', 100000);

mbsys = mbdyn.pre.system ( {prb}, ...
                           'Nodes', nodes, ...
                           'Elements', [elements, bodies, {pinjoint, omega_joint, mbdyn.pre.gravity()} ], ...
                           'DefaultOrientation', 'euler123');
                       
str = mbsys.generateMBDynInputStr ()

% mbsys.setStructuralNodeSize (L/10, L/10, L/10);
% 
% mbsys.draw ( 'Mode', 'wireghost', ...
%              'References', true, ...
%              'ReferenceScale', 0.5, ...
%              'AxLims', [-2, 2; -2, 2; -2, 2;] );
         
%% Generate input file and run mbdyn

inputfile = mbsys.generateMBDynInputFile ('Test_mbdyn_axialRotation.mbd');

% create the mbdyn input file
mbsys.generateMBDynInputFile (inputfile);

% start mbdyn 
mbdyn.mint.start_mbdyn (inputfile);

%% Post-processing

% load data into post-processing object
mbdynpost = mbdyn.postproc (inputfile(1:end-4), mbsys);

%% Animate

% mbdynpost.animate ( 'PlotTrajectories', false, ...
%                     'DrawLabels', false, ...
%                     ...'Skip', 20, ...
%                     'DrawMode', 'wireghost', ...
%                     'Light', true, ...
%                     'AxLims', [-2, 2; -2, 2; -2, 2;], ...
%                     'VideoFile', 'Test_mbdyn_magnetic_gear.avi');
                
%% Plot angular velocities

mbdynpost.plotNodeAngularVelocities ();
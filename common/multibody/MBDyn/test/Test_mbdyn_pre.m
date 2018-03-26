%% mbdyn.pre.matrixFreeSolver('bicgstab')

slv = mbdyn.pre.matrixFreeSolver('bicgstab');

str = slv.generateMBDynInputString ()

%% linearsolver

lslv = mbdyn.pre.linearSolver ('umfpack');

str = lslv.generateMBDynInputString ()

lslv = mbdyn.pre.linearSolver ('umfpack', 'BlockSize', 32);

str = lslv.generateMBDynInputString ()

%% initial value

itime = 0;
ftime = 1; 
tstep = 0.1;
pbm = mbdyn.pre.initialValueProblem (itime, ftime, tstep, 'Output', {'iterations', 'residual'});

str = pbm.generateMBDynInputString ()

pbm = mbdyn.pre.initialValueProblem (itime, ftime, tstep, ...
                    'Output', {'iterations', 'residual'}, ...
                    'ResidualTolerance', 1e-6);

str = pbm.generateMBDynInputString ()

pbm = mbdyn.pre.initialValueProblem (itime, ftime, tstep, ...
                    'Output', {'iterations', 'residual'}, ...
                    'ResidualTolerance', {1e-6, 'test', 'minmax', 'scale'}, ...
                    'SolutionTolerance', {1e-6, 'test', 'minmax'});

str = pbm.generateMBDynInputString ()

%% orientation matrix

om = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1;0;0], 'ib', 2, 'vecB', [0;1;0])); 

om.orientationMatrix


%% orientation matrix

pos = [1,0,0]

om = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [cos(pi/6);sin(pi/6);0], 'ib', 3, 'vecB', [0;0;1])); 

om.orientationMatrix

pos * om.orientationMatrix

om = mbdyn.pre.orientmat ('euler', [0,0,pi/6]);

om.orientationMatrix

pos * om.orientationMatrix

% plot the orientation matrix as 3 axes
om.draw ();

om = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1,0,0], 'ib', 2, 'vecB', 'guess')); 

om.draw ()

om = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [0.5,0.5,0], 'ib', 2, 'vecB', [])); 

om.draw ()

%% references

gref = mbdyn.pre.globalref

theta1 = 0.0;
theta2 = pi/4;
L = 1;

Ref_Link1 = mbdyn.pre.reference ([], mbdyn.pre.orientmat ('euler', [0, pi/2 - theta1, 0]), [], [], 'Parent', gref)

Ref_Link1.pos

Ref_Link2 = mbdyn.pre.reference ([L; 0; 0], mbdyn.pre.orientmat ('euler', [0, -theta2, 0]), [], [], 'Parent', Ref_Link1)

Ref_Link2.pos

Node_Link2 = mbdyn.pre.reference ([0.5*L; 0; 0], mbdyn.pre.orientmat ('orientation', eye(3)), [], [], 'Parent', Ref_Link2)

fcalc = [ L*sin(theta1)+(L/2)*sin(theta1+theta2); 
          0; 
          -L*cos(theta1)-(L/2)*cos(theta1+theta2) ];

      
[ fcalc, Node_Link2.pos ]

% plot thre references
mbdyn.pre.base.drawReferences ({Ref_Link1, Ref_Link2}, 'Scale', 0.5)

% om = mbdyn.pre.orientmat ('euler', [0, pi-(theta1+theta2), 0])

%% references convertGlobal

ref1 = mbdyn.pre.reference ([0;1;0], mbdyn.pre.orientmat ('euler123', [0,pi/4,pi/4]), [1;0;0], []);

pos = [1;1;0];
orientm = mbdyn.pre.orientmat ('euler123', [0,0,pi/2]);
vel = [1;0;0];
omega = [0;0.5;0];

ref2 = mbdyn.pre.reference (pos, orientm, [], []);

[dpos, dorientm, dvel, domega] = ref1.convertGlobal (pos, orientm, vel, omega)

ref3 = mbdyn.pre.reference (dpos, dorientm, dvel, domega, ...
                            'Parent', ref1);

mbdyn.pre.base.drawReferences ( {ref1, ref2, ref3} );

[ pos, ref3.pos ]
[ vel, ref3.v ]
orientm.orientationMatrix
ref3.orientm.orientationMatrix
[ omega, ref3.omega ]

%% base node

abn = mbdyn.pre.node ();

str = abn.generateMBDynInputString ()


abn = mbdyn.pre.node ('Scale', 3);

str = abn.generateMBDynInputString ()

abn = mbdyn.pre.node ('Scale', 'default');

str = abn.generateMBDynInputString ()

try
abn = mbdyn.pre.node ('Scale', 'fdasfsa');
catch
    disp ('Correctly caught bad Scale string');
end

abn = mbdyn.pre.node ('Output', 'no');

str = abn.generateMBDynInputString ()

abn = mbdyn.pre.node ('Scale', 3, 'Output', 'yes');

str = abn.generateMBDynInputString ()

try
abn = mbdyn.pre.node ('Output', 'fdasfsa');
catch
    disp ('Correctly caught bad Output string');
end

try
abn = mbdyn.pre.node ('Output', 9);
catch
    disp ('Correctly caught bad Output value');
end

%% abstract node

abn = mbdyn.pre.abstractNode ();

str = abn.generateMBDynInputString ()

abn = mbdyn.pre.abstractNode ('Value', 10);

str = abn.generateMBDynInputString ()

abn = mbdyn.pre.abstractNode  ('Value', 10, 'Derivative', 100);

str = abn.generateMBDynInputString ()

abn = mbdyn.pre.abstractNode  ('Derivative', 100)

str = abn.generateMBDynInputString ()

abn = mbdyn.pre.abstractNode ('Value', 10, 'Scale', 3, 'Output', 'yes');

str = abn.generateMBDynInputString ()

%% structuralNode6dof


sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

str = sn6dof.generateMBDynInputString ()

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true, ...
                        'Scale', 3, 'Output', 'yes');
                    
str = sn6dof.generateMBDynInputString ()

%% structuralNode3dof

sn3dof = mbdyn.pre.structuralNode3dof ('dynamic displacement', 'Accel', true);

str = sn3dof.generateMBDynInputString ()

%% structuralNodeDummy

sn3dof = mbdyn.pre.structuralNode6dof ('dynamic displacement', 'Accel', true);

str = sn3dof.generateMBDynInputString ()

%% Body

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
mass = 1;
cog = [0;0;0];
inertiamat = eye (3);

bd = mbdyn.pre.body (mass, cog, inertiamat, sn6dof);

str = bd.generateMBDynInputString ()

%% Body

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
mass = 1;
cog = [0;0;0];
inertiamat = eye (3);

bd = mbdyn.pre.body (mass, cog, inertiamat, sn6dof, 'InertialOrientation', eye (3));

str = bd.generateMBDynInputString ()

%% Body

sn3dof = mbdyn.pre.structuralNode3dof ('dynamic displacement', 'Accel', true);

mass = 1;

bd = mbdyn.pre.body (mass, [], [], sn3dof);

str = bd.generateMBDynInputString ()

%% Multiple Mass Body

% body: BODY_LABEL, NODE_LABEL,
% condense, 3,
%   4., # mass 1 (mid)
%   reference, node, 0., 0., 0., # c.m. offset 1
%   diag, .4, .4, .2, # inertia tensor 1
%   2., # mass 2 (top)
%   reference, node, 0., 0., 1., # c.m. offset 2
%   diag, .2, .2, .1, # inertia tensor 2
%   2., # mass 3 (bottom)
%   reference, node, 0., 0., -1., # c.m. offset 3
%   diag, .2, .2, .1; # inertia tensor 3

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
mass = [4, 2, 2];
cog = { [0;0;0], [0;0;1], [0;0;-1] };
inertiamat = { diag([0.4, 0.4, 0.2]), diag([0.2,0.2,0.1]), diag([0.2,0.2,0.1]) };

bd = mbdyn.pre.bodyMultiMass (mass, cog, inertiamat, sn6dof);

str = bd.generateMBDynInputString ()

bd.setSize (1, 0.2, 0.2, 0.2);
bd.setSize (2, 0.1, 0.1, 0.1);
bd.setSize (3, 0.1, 0.1, 0.1);
bd.draw ()


%% Total Joint

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

posstatus = 'active';
orientstatus = true;

jnt = mbdyn.pre.totalJoint (sn1, sn2, posstatus, orientstatus);
jnt.generateMBDynInputString ()


jnt = mbdyn.pre.totalJoint (sn1, sn2, posstatus, orientstatus, ...
    'RelativeOffset1', [1; 2; 3]);
jnt.generateMBDynInputString ()

jnt = mbdyn.pre.totalJoint (sn1, sn2, posstatus, orientstatus, ...
    'RelativeOffset1', [1; 2; 3], ...
    'RelativeOffset1Reference', 'other node');
jnt.generateMBDynInputString ()

%% Revolute Rotation

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

jnt = mbdyn.pre.revoluteRotation (sn1, sn2);
jnt.generateMBDynInputString ()

%% Revolute Hinge

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

friction_fcn = mbdyn.pre.constScalarFunction ('constant roller bearing', 0.0018);

friction_model = mbdyn.pre.discreteCoulombFriction (friction_fcn);
friction_radius = 1;
radius = friction_radius;
shp = mbdyn.pre.simpleShape ();

revhinge = mbdyn.pre.revoluteHinge ( sn1, sn2, ...
                                     [0;0;0], [0;0;0], ...
                                     'Offset1Reference', 'node', ...
                                     'Offset2Reference', 'other node', ...
                                     'RelativeOrientation1', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation1Reference', 'node', ...
                                     'RelativeOrientation2', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation2Reference', 'other node', ...
                                     'FrictionRadius', friction_radius, ...
                                     'FrictionModel', friction_model ...
                                     , 'ShapeFUnction', shp ...
                                   );
                               
revhinge.generateMBDynInputString ()

%%
preload = 1;
revhinge = mbdyn.pre.revoluteHinge ( sn1, sn2, ...
                                     [0;0;0], [0;0;0], ...
                                     'Offset1Reference', 'node', ...
                                     'Offset2Reference', 'other node', ...
                                     'RelativeOrientation1', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation1Reference', 'node', ...
                                     'RelativeOrientation2', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation2Reference', 'other node', ...
                                     'FrictionRadius', friction_radius, ...
                                     'FrictionModel', friction_model ...
                                     , 'Preload', preload ...
                                   );
                               

                               
revhinge.generateMBDynInputString ()
                                 
%% system

% this is intended to replicate the sky-engineering example at:
%
% http://www.sky-engin.jp/en/MBDynTutorial/chap15/chap15.html

gref = mbdyn.pre.globalref;

theta1 = -pi/5;
theta2 = pi/5;
L = 1;
M = 1;
inertiamat = diag ([0., M*L^2./12., M*L^2./12.]);

Ref_Link1 = mbdyn.pre.reference ( [], ...
                                  mbdyn.pre.orientmat ('euler', [0, theta1, 0]), ...
                                  [], ...
                                  [], ...
                                  'Parent', gref, ...
                                  'Name', 'Ref_Link1');
                              

Ref_Node_Link1 = mbdyn.pre.reference ( [0.5*L; 0; 0], ...
                                       [], ...
                                       [], ...
                                       [], ...
                                       'Parent', Ref_Link1, ...
                                       'Name', 'Ref_Node_Link1' );

link1node = mbdyn.pre.structuralNode6dof ('dynamic', ...
                                          'AbsolutePosition', Ref_Node_Link1.pos, ...
                                          'AbsoluteOrientation', Ref_Node_Link1.orientm );

Ref_Link2 = mbdyn.pre.reference ( [L; 0; 0], ...
                                   mbdyn.pre.orientmat ('euler', [0, theta2, 0]), [], ...
                                   [], ...
                                   'Parent', Ref_Link1, ...
                                   'Name', 'Ref_Link2');

Ref_Node_Link2 = mbdyn.pre.reference ( [0.5*L; 0; 0], ...
                                       mbdyn.pre.orientmat ('orientation', eye(3)), ...
                                       [], ...
                                       [], ...
                                       'Parent', Ref_Link2, ...
                                       'Name', 'Ref_Node_Link2');

link2node = mbdyn.pre.structuralNode6dof ('dynamic', ...
                                          'AbsolutePosition', Ref_Node_Link2.pos, ...
                                          'AbsoluteOrientation', Ref_Node_Link2.orientm );

link1 = mbdyn.pre.body (M, [], inertiamat, link1node);
link2 = mbdyn.pre.body (M, [], inertiamat, link2node);

link1.setSize (L, L/10, L/10);
link2.setSize (L, L/10, L/10);
link1.setColour ('r');
link1.setColour ('b');

hinges_orientation = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1;0;0], 'ib', 3, 'vecB', [0;1;0]));

Ref_pin = mbdyn.pre.reference ([], hinges_orientation, [], [], 'Parent', Ref_Link1, 'Name', 'Ref_pin');
Ref_hinge = mbdyn.pre.reference ([], hinges_orientation, [], [], 'Parent', Ref_Link2, 'Name', 'Ref_hinge');

pinjoint = mbdyn.pre.revolutePin (link1node, Ref_Link1.pos, Ref_Link1.pos, ...
                    'NodeRelativeOrientation', Ref_pin.orientm, ...
                    'PinOrientation', Ref_pin.orientm);
                
linkjoint = mbdyn.pre.revoluteHinge ( link1node, link2node, ...
                                      ...Ref_Link2.pos, Ref_Link2.pos, ...
                                     [L/2;0;0], [L/2;0;0], ...
                                     'Offset1Reference', 'node', ...
                                     'Offset2Reference', 'other node', ...
                                     ...'Offset1Reference', 'global', ...
                                     ...'Offset2Reference', 'global', ...
                                     'RelativeOrientation1', Ref_hinge.orientm, ...
                                     'Orientation1Reference', 'global', ...
                                     'RelativeOrientation2', Ref_hinge.orientm, ...
                                     'Orientation2Reference', 'global' ...
                                     );
                
pinjoint.setSize (L/10, L/10, 2*L/10);
pinjoint.setColour ('k');
linkjoint.setSize (L/10, L/10, 2*L/10);
linkjoint.setColour ('g');

prb = mbdyn.pre.initialValueProblem (0, 10, 1e-3, 'ResidualTolerance', 1e-9) %, 'DerivativesTolerance', 100000);

mbsys = mbdyn.pre.system ( {prb}, ...
                           'Nodes', {link1node, link2node}, ...
                           'Elements', {link1, link2, pinjoint, linkjoint, mbdyn.pre.gravity()}, ...
                           'References', {Ref_pin, Ref_hinge}, ... {Ref_Link1, Ref_Node_Link1, Ref_Link2, Ref_Node_Link2, Ref_pin, Ref_hinge}
                           'DefaultOrientation', 'euler123');

str = mbsys.generateMBDynInputStr ()

mbsys.setStructuralNodeSize (L/10, L/10, L/10);

mbsys.draw ( 'Mode', 'wireframe', ...
             'References', false, ...
             'ReferenceScale', 0.5, ...
             'AxLims', [-2, 2; -2, 2; -2, 2;] );

%% Generate input file and run mbdyn

inputfile = mbsys.generateMBDynInputFile ('Test_mbdyn_pre.mbd');

% create the mbdyn input file
mbsys.generateMBDynInputFile (inputfile);

% start mbdyn 
mbdyn.mint.start_mbdyn (inputfile);

%% Post-processing

% load data into post-processing object
mbdynpost = mbdyn.postproc (inputfile(1:end-4), mbsys);

%% Plot Trajectories

mbdynpost.plotNodeTrajectories ();

%% Plot a particular time step of interest

mbdynpost.drawStep (500, 'AxLims', [-2, 2; -2, 2; -2, 2;]);

%% Animate

mbdynpost.animate ( 'PlotTrajectories', true, ...
                    'DrawLabels', true, ...
                    'Skip', 20, ...
                    'DrawMode', 'solid', ...
                    'Light', true, ...
                    'AxLims', [-3.1, 3.1; -1.5, 1.5;  -3, 3], ...
                    'VideoFile', 'double_pendulum.avi');

%% Drawing

el = mbdyn.pre.element ();
% 
el.draw ()
el.draw('Mode', 'wireframe');

%%

% om = mbdyn.pre.orientmat ('2vectors', struct ('vec1axis', 1, 'vec1', [cos(pi/6);sin(pi/6);0], 'vec2axis', 3, 'vec2', [0;0;1])); 

om = mbdyn.pre.orientmat ('euler', [0,0,pi/6]);

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', ...
                                       'AbsolutePosition', [0.5; 0.5; 0], ...
                                       'AbsoluteOrientation', om, ...
                                       'Accel', true);
mass = 1;
cog = [0;0;0];
inertiamat = eye (3);

bd = mbdyn.pre.body (mass, cog, inertiamat, sn6dof);
bd.draw ('Mode', 'wireghost')
xlabel ('x'); ylabel ('y'); zlabel('z'); view (3)

%% socket communicator

soc = mbdyn.pre.socketCommunicator ('Path', '/tmp/mbdyn.sock');
soc.generateMBDynInputString ()

%% external sturctural force

soc = mbdyn.pre.socketCommunicator ('Path', '/tmp/mbdyn.sock');
soc.generateMBDynInputString ()


sn6dof1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn6dof2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true, 'AbsolutePosition', [1;0;0]);
sn6dof3 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true, 'AbsolutePosition', [2;0;0]);

nodeoffsets = struct ('NodeInd', {1, 3}, ...
                      'Offset', {[0;1;0], [1;1;1]}, ...
                      'OffsetType', {'global', 'local'})
                  
extsf = mbdyn.pre.externalStructuralForce ( {sn6dof1, sn6dof2, sn6dof3}, ...
                                            nodeoffsets, ...
                                            soc);
                                        
extsf.generateMBDynInputString ()



%% stateSpaceFilter

state_order = 2;
A = [ 1, 2; 3, 4 ];
B = [ 1, 2; 3, 4 ] * 2;
C = [ 1, 2; 3, 4 ] * 3;
D = [ 1, 2; 3, 4 ] * 4;
E = [ 1, 2; 3, 4 ] * 5;

SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

%%
SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'gain', 1, ...
                                  'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'D', D);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'E', E);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'D', D, 'E', E);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');



SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'D', D, 'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'E', E, 'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'D', D, 'E', E, 'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');




SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'D', D, 'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'E', E, 'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceFilter (state_order, A, B, C, ...
                                  'D', D, 'E', E, 'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


%% nodeDOF

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

dof = mbdyn.pre.nodeDOF(sn1, 'DOFNumber', 2, 'Alge', 'algebraic');

dof.generateMBDynInputString ()

sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
dof2 = mbdyn.pre.nodeDOF(sn1, 'DOFNumber', 1, 'Alge', 'algebraic');

isa ([dof, dof2], 'mbdyn.pre.nodeDOF')

%% stateSpaceMIMO

state_order = 2;
A = [ 1, 2; 3, 4 ];
B = [ 1, 2; 3, 4 ] * 2;
C = [ 1, 2; 3, 4 ] * 3;
D = [ 1, 2; 3, 4 ] * 4;
E = [ 1, 2; 3, 4 ] * 5;

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

dof1 = mbdyn.pre.nodeDOF(sn1, 'DOFNumber', 2, 'Alge', 'algebraic');

sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
dof2 = mbdyn.pre.nodeDOF(sn1, 'DOFNumber', 1, 'Alge', 'algebraic');

output_node_list = [dof1, dof2];
input = {dof1, dof2};

SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'gain', 1, ...
                                  'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'D', D);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'E', E);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'D', D, 'E', E);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');



SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'D', D, 'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'E', E, 'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'D', D, 'E', E, 'gain', 1);

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');




SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'D', D, 'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'E', E, 'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');

SSF = mbdyn.pre.stateSpaceMIMO (state_order, A, B, C, output_node_list, input, ...
                                  'D', D, 'E', E, 'balance', 'no');

SSF.generateMBDynInputString ()
fprintf (1, '\n\n');


%% nodeDrive

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

func_drive = mbdyn.pre.directDrive ();

nd = mbdyn.pre.nodeDrive (sn1, func_drive, 'Index', 2);

nd.generateMBDynInputString ()

nd = mbdyn.pre.nodeDrive (sn1, func_drive,  'String', 'XP[1]');

nd.generateMBDynInputString ()


%% compnent template drive caller

drivecallers = {mbdyn.pre.const(1), 'inactive', mbdyn.pre.const(2)};

dc = mbdyn.pre.componentTplDriveCaller (drivecallers);

dc.generateMBDynInputString ()

dc = mbdyn.pre.componentTplDriveCaller (drivecallers, 'ShapeType', 'sym');

dc.generateMBDynInputString ()

dc = mbdyn.pre.componentTplDriveCaller (drivecallers, 'ShapeType', 'diag');

dc.generateMBDynInputString ()

dc = mbdyn.pre.componentTplDriveCaller (drivecallers, 'ShapeType', 'bum');

dc.generateMBDynInputString ()


%% total force

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
drivecallers = {mbdyn.pre.const(1), 'inactive', mbdyn.pre.const(2)};
dc = mbdyn.pre.componentTplDriveCaller (drivecallers);

tf = mbdyn.pre.totalForce  (sn1, 'Force', dc);

tf.generateMBDynInputString ()

tf = mbdyn.pre.totalForce  (sn1, 'Force', dc, 'Moment', dc);

tf.generateMBDynInputString ()


newabsnodes = {mbdyn.pre.abstractNode(), mbdyn.pre.abstractNode(), mbdyn.pre.abstractNode()};

drivecallers = { mbdyn.pre.nodeDrive(newabsnodes{1}, mbdyn.pre.directDrive()), ...
                 mbdyn.pre.nodeDrive(newabsnodes{2}, mbdyn.pre.directDrive()), ...
                 mbdyn.pre.nodeDrive(newabsnodes{3}, mbdyn.pre.directDrive()) };

dc = mbdyn.pre.componentTplDriveCaller (drivecallers);

tf = mbdyn.pre.totalForce  (sn1, 'Moment', dc);

tf.generateMBDynInputString ()

%% simpleShape

shp = mbdyn.pre.simpleShape ();

shp.generateMBDynInputString ()

%% simplePlaneHingeShape

radius = pi;
shp = mbdyn.pre.simplePlaneHingeShape (radius);

shp.generateMBDynInputString ()

%% const scalar function

x = mbdyn.pre.constScalarFunction ('jimbob', tau);

x.generateMBDynInputString ()

%% discreteCoulombFriction Model

friction_fcn = mbdyn.pre.constScalarFunction ('jimbob', tau);

x = mbdyn.pre.discreteCoulombFriction (friction_fcn);

x.generateMBDynInputString ()


%% angularVelocity

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
w = 1;
drv = mbdyn.pre.const (w);
rotaxis = [1;0;0];

av = mbdyn.pre.angularVelocity ( sn1, rotaxis, drv );

av.generateMBDynInputString ()


%% rampDrive

rd = mbdyn.pre.rampDrive (1, 2, 4, 10);

rd.generateMBDynInputString ()

rd = mbdyn.pre.rampDrive (1, 'forever', 4, 10);

rd.generateMBDynInputString ()

%% stringDrive

sd = mbdyn.pre.stringDrive ('5 * Var');

sd.generateMBDynInputString ()

%% deformableAxialJoint

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

def_fcn = mbdyn.pre.sinScalarFunction ('sin wave', 2, 1);

law = mbdyn.pre.scalarFunctionElasticIsotropicConstituativeLaw (def_fcn);

obj = mbdyn.pre.deformableAxialJoint ( sn1, sn2, law, ...
                                     'Offset1', [0;0;0], ...
                                     'Offset2', [0;0;0], ...
                                     'Offset1Reference', 'node', ...
                                     'Offset2Reference', 'other node', ...
                                     'RelativeOrientation1', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation1Reference', 'node', ...
                                     'RelativeOrientation2', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation2Reference', 'other node' ...
                                   );

obj.generateMBDynInputString ()


%% axialRotation

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

t_init = 0;
t_final = 1;
slope = 1;
init_val = 0; 

omega_drive = mbdyn.pre.rampDrive (t_init, t_final, slope, init_val);

obj = mbdyn.pre.axialRotation ( sn1, sn2, omega_drive, ...
                                     'Offset1', [0;0;0], ...
                                     'Offset2', [0;0;0], ...
                                     'Offset1Reference', 'node', ...
                                     'Offset2Reference', 'other node', ...
                                     'RelativeOrientation1', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation1Reference', 'node', ...
                                     'RelativeOrientation2', mbdyn.pre.orientmat ('eye'), ...
                                     'Orientation2Reference', 'other node' ...
                                   );

obj.generateMBDynInputString ()

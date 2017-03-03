

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

%% references

gref = mbdyn.pre.globalref

theta1 = 0.1;
theta2 = 0.1;
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


% om = mbdyn.pre.orientmat ('euler', [0, pi-(theta1+theta2), 0])

%% structuralNode6dof


sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

str = sn6dof.generateOutputString ()

%%
sn3dof = mbdyn.pre.structuralNode3dof ('dynamic displacement', 'Accel', true);

str = sn3dof.generateOutputString ()

%% Body

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
mass = 1;
cog = [0;0;0];
inertiamat = eye (3);

bd = mbdyn.pre.body (mass, cog, inertiamat, sn6dof);

str = bd.generateOutputString ()

%%

sn6dof = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
mass = 1;
cog = [0;0;0];
inertiamat = eye (3);

bd = mbdyn.pre.body (mass, cog, inertiamat, sn6dof, 'InertialOrientation', eye (3));

str = bd.generateOutputString ()


%% Total Joint

sn1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);

posstatus = 'active';
orientstatus = true;

jnt = mbdyn.pre.totalJoint (sn1, sn2, posstatus, orientstatus);
jnt.generateOutputString ()


jnt = mbdyn.pre.totalJoint (sn1, sn2, posstatus, orientstatus, ...
    'RelativeOffset1', [1; 2; 3]);
jnt.generateOutputString ()

jnt = mbdyn.pre.totalJoint (sn1, sn2, posstatus, orientstatus, ...
    'RelativeOffset1', {'other node', [1; 2; 3]});
jnt.generateOutputString ()


%% system

% this is intended to replicate the sky-engineering example at:
%
% http://www.sky-engin.jp/en/MBDynTutorial/chap15/chap15.html

gref = mbdyn.pre.globalref;

theta1 = pi/2;
theta2 = 0;
L = 1;
M = 1;
inertiamat = diag ([0., M*L^2./12., M*L^2./12.]);

Ref_Link1 = mbdyn.pre.reference ( [], ...
                                  mbdyn.pre.orientmat ('euler', [0, pi/2 - theta1, 0]), ...
                                  [], ...
                                  [], ...
                                  'Parent', gref);

Ref_Node_Link1 = mbdyn.pre.reference ( [0.5*L; 0; 0], ...
                                       [], ...
                                       [], ...
                                       [], ...
                                       'Parent', Ref_Link1 );

link1node = mbdyn.pre.structuralNode6dof ('dynamic', ...
                                          'AbsolutePosition', Ref_Node_Link1.pos, ...
                                          'AbsoluteOrientation', Ref_Node_Link1.orientm );

Ref_Link2 = mbdyn.pre.reference ( [L; 0; 0], ...
                                   mbdyn.pre.orientmat ('euler', [0, -theta2, 0]), [], ...
                                   [], ...
                                   'Parent', Ref_Link1 );

Ref_Node_Link2 = mbdyn.pre.reference ( [0.5*L; 0; 0], ...
                                       mbdyn.pre.orientmat ('orientation', eye(3)), ...
                                       [], ...
                                       [], ...
                                       'Parent', Ref_Link2);

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

Ref_pin = mbdyn.pre.reference ([], hinges_orientation, [], [], 'Parent', Ref_Link1);
Ref_hinge = mbdyn.pre.reference ([], hinges_orientation, [], [], 'Parent', Ref_Link2);

pinjoint = mbdyn.pre.revolutePin (link1node, Ref_Link1.pos, Ref_Link1.pos, ...
                    'NodeRelativeOrientation', hinges_orientation, ...
                    'PinOrientation', hinges_orientation);
                
linkjoint = mbdyn.pre.revoluteHinge (link1node, link2node, Ref_Link2.pos, Ref_Link2.pos, ...
                    'Offset1Reference', 'global', ...
                    'Offset2Reference', 'global', ...
                    'RelativeOrientation1', Ref_hinge.orientm, ...
                    'Orientation1Reference', 'global', ...
                    'RelativeOrientation2', Ref_hinge.orientm, ...
                    'Orientation2Reference', 'global');
                
pinjoint.setSize (L/10, L/10, L/10);
pinjoint.setColour ('k');
linkjoint.setSize (L/10, L/10, L/10);
linkjoint.setColour ('g');

prb = mbdyn.pre.initialValueProblem (0, 5, 1e-3);

mbsys = mbdyn.pre.system ( {prb}, ...
                           'Nodes', {link1node, link2node}, ...
                           'Elements', {link1, link2, pinjoint, linkjoint, mbdyn.pre.gravity()} );

str = mbsys.generateMBDynInputStr ()

mbsys.setStructuralNodeSize (L/10, L/10, L/10);
mbsys.draw ('Mode','wireghost')


filename = mbsys.generateMBDynInputFile ('Test_mbdyn_pre.mbd');

% start mbdyn 
% delete ('output.*');

[status, cmdout] = system ( sprintf ('export LD_LIBRARY_PATH="" ; mbdyn -f "%s" -o "%s" > "%s" 2>&1 &', ...
                    filename, ...
                    [filename(1:end-4), '_mbd'], ...
                    [filename(1:end-4), '_mbd.txt'] ...
                                     ) ...
                           );
                     
%% Post-processing

mbdynpost = mbdyn.postproc (mbsys);
mbdynpost.loadResultsFromFiles ([filename(1:end-4), '_mbd']);

mbdynpost.plotNodeTrajectories ()

mbdynpost.animate ('PlotTrajectories', true, 'DrawLabels', true, 'Skip', 40)

% [status, cmdout] = system (sprintf ('mbdyn -f "%s" -o output > output.txt 2>&1 &', mbdpath))

% linkjoint.draw ()

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
soc.generateOutputString ()

%% external sturctural force

soc = mbdyn.pre.socketCommunicator ('Path', '/tmp/mbdyn.sock');
soc.generateOutputString ()


sn6dof1 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true);
sn6dof2 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true, 'AbsolutePosition', [1;0;0]);
sn6dof3 = mbdyn.pre.structuralNode6dof ('dynamic', 'Accel', true, 'AbsolutePosition', [2;0;0]);

nodeoffsets = struct ('NodeInd', {1, 3}, ...
                      'Offset', {[0;1;0], [1;1;1]}, ...
                      'OffsetType', {'global', 'local'})
                  
extsf = mbdyn.pre.externalStructuralForce ( {sn6dof1, sn6dof2, sn6dof3}, ...
                                            nodeoffsets, ...
                                            soc);
                                        
extsf.generateOutputString ()
                                        
                                        
                                        






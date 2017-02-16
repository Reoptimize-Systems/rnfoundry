

%% orientation matrix

om = mbdyn.pre.orientmat ('2vectors', struct ('vec1axis', 1, 'vec1', [1;0;0], 'vec2axis', 2, 'vec2', [0;1;0])); 

om.orientationMatrix


%% orientation matrix
om = mbdyn.pre.orientmat ('2vectors', struct ('vec1axis', 1, 'vec1', [cos(pi/6);sin(pi/6);0], 'vec2axis', 3, 'vec2', [0;0;1])); 

om.orientationMatrix

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

%%

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

gref = mbdyn.pre.globalref;

theta1 = 0.1;
theta2 = 0.1;
L = 1;
M = 1;
inertiamat = diag ([0., M*L^2./12., M*L^2./12.]);

Ref_Link1 = mbdyn.pre.reference ([], mbdyn.pre.orientmat ('euler', [0, pi/2 - theta1, 0]), [], [], 'Parent', gref);

Ref_Node_Link1 = mbdyn.pre.reference ([0.5*L; 0; 0], [], [], [], 'Parent', Ref_Link1);

link1node = mbdyn.pre.structuralNode6dof ('dynamic', 'AbsolutePosition', Ref_Link1.pos);

Ref_Link2 = mbdyn.pre.reference ([L; 0; 0], mbdyn.pre.orientmat ('euler', [0, -theta2, 0]), [], [], 'Parent', Ref_Link1);

Ref_Node_Link2 = mbdyn.pre.reference ([0.5*L; 0; 0], mbdyn.pre.orientmat ('orientation', eye(3)), [], [], 'Parent', Ref_Link2);

link2node = mbdyn.pre.structuralNode6dof ('dynamic', 'AbsolutePosition', Ref_Link2.pos);


link1 = mbdyn.pre.body (M, [], inertiamat, link1node);
link2 = mbdyn.pre.body (M, [], inertiamat, link2node);

hinges_orientation = mbdyn.pre.orientmat ('2vectors', struct ('vec1axis', 1, 'vec1', [1;0;0], 'vec2axis', 3, 'vec2', [0;1;0]));

Ref_pin = mbdyn.pre.reference ([], hinges_orientation, [], [], 'Parent', Ref_Link1);
Ref_hinge = mbdyn.pre.reference ([], hinges_orientation, [], [], 'Parent', Ref_Link2);

pinjoint = mbdyn.pre.revolutePin (link1node, Ref_Link1.pos, Ref_Link1.pos, ...
                    'NodeRelativeOrientation', hinges_orientation, ...
                    'PinOrientation', hinges_orientation);
                
linkjoint = mbdyn.pre.revoluteHinge (link1node, link2node, Ref_Link2.pos, Ref_Link2.pos, ...
                    'RelativeOrientation1', Ref_hinge.orientm, ...
                    'RelativeOrientation2', Ref_hinge.orientm );

mbsys = mbdyn.pre.system ( 'Nodes', {link1node, link2node}, ...
                           'Elements', {link1, link2, pinjoint, linkjoint} );

str = mbsys.generateMBDynInputStr ()






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
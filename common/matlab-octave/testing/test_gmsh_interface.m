
bobj = gmsh.base ();

bobj.geoName = 'Test';
bobj.tag = 1;

bobj.generateGeoFileStr ()

%% point

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0, 'ElementSize', 0.1);

pnt1.tag = 1;
pnt2.tag = 2;

pnt1.generateGeoFileStr ()

pnt2.generateGeoFileStr ()


%% physical point

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0);

pnt1.tag = 1;
pnt2.tag = 2;


ppnt1 = gmsh.physicalPoint ([pnt1, pnt2]);

ppnt1.tag = 3;

ppnt1.generateGeoFileStr ()


ppnt2 = gmsh.physicalPoint ([pnt1, pnt2], 'Label', 'test label');

ppnt2.tag = 4;

ppnt2.generateGeoFileStr ()


%% line

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0);

pnt1.tag = 1;
pnt2.tag = 2;


line = gmsh.line (pnt1, pnt2);

line.tag = 3;

line.generateGeoFileStr ()

%% bezier

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0);

pnt3 = gmsh.point (2,0,0);

pnt4 = gmsh.point (3,0,0);

pnt1.tag = 1;
pnt2.tag = 2;
pnt3.tag = 3;
pnt4.tag = 4;

bez = gmsh.bezier ([pnt1, pnt2, pnt3, pnt4]);

bez.tag = 5;

bez.generateGeoFileStr ()

%% bspline

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0);

pnt3 = gmsh.point (2,0,0);

pnt4 = gmsh.point (3,0,0);

pnt1.tag = 1;
pnt2.tag = 2;
pnt3.tag = 3;
pnt4.tag = 4;

bspln = gmsh.bSpline ([pnt1, pnt2, pnt3, pnt4]);

bspln.tag = 5;

bspln.generateGeoFileStr ()

%% spline

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0);

pnt3 = gmsh.point (2,0,0);

pnt4 = gmsh.point (3,0,0);

pnt1.tag = 1;
pnt2.tag = 2;
pnt3.tag = 3;
pnt4.tag = 4;

spln = gmsh.spline ([pnt1, pnt2, pnt3, pnt4]);

spln.tag = 5;

spln.generateGeoFileStr ()

%% circle

pnt1 = gmsh.point (0,0,0);

pnt2 = gmsh.point (1,0,0);

pnt3 = gmsh.point (0.5,0,0);

pnt1.tag = 1;
pnt2.tag = 2;
pnt3.tag = 3;

circle = gmsh.circle (pnt1, pnt2, pnt3);

circle.tag = 4;

circle.generateGeoFileStr ()


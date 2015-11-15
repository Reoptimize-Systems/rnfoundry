% Test_axialfluxstatorhalf2dfemmprob


slots = 2 * 21;
Poles = 28;
ypole = 1;
yslot = ypole / (slots / Poles);
ycoil = yslot * 0.8;
yshoegap = 0.2 * ycoil;
xyoke = 0.05;
xcoil = 0.3;
xshoebase = 0.05;
xshoegap = 0.025;
xoffset = 0;
side = 'l';
nslots = 3;

[FemmProblem, outernodes, coillabellocs, slotinfo] = ...
    uncurvedstatorhalf2dfemmprob(nslots, yslot, ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, xoffset, side, 'NWindingLayers', 3);

plotfemmproblem (FemmProblem);

%%
filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))


%% test any number of slots

slots = 2 * 21;
Poles = 28;
ypole = 1;
yslot = ypole / (slots / Poles);
ycoil = yslot * 0.8;
yshoegap = 0.2 * ycoil;
xyoke = 0.05;
xcoil = 0.3;
xshoebase = 0.05;
xshoegap = 0.025;
xoffset = 0;
side = 'l';

[FemmProblem, outernodes, coillabellocs] = ...
    uncurvedstatorhalf2dfemmprob(nslots, yslot, ypole, ycoil, yshoegap, ...
                                  xyoke, xcoil, xshoebase, xshoegap, xoffset, side, ...
                                  'NWindingLayers', 2, ...
                                  'NSlots', 4);

filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))
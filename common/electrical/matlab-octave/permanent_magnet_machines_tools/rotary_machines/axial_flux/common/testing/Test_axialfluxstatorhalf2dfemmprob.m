% Test_axialfluxstatorhalf2dfemmprob


slots = 2 * 21;
poles = 28;
ypole = 1;
ycoil = ypole / (slots / poles) * 0.8;
yshoegap = 0.2 * ycoil;
xyoke = 0.05;
xcoil = 0.3;
xshoebase = 0.05;
xshoegap = 0.025;
xoffset = 0;
side = 'l';

[FemmProblem, outernodes] = ...
    axialfluxstatorhalf2dfemmprob(slots, poles, ypole, ycoil, yshoegap, ...
                                  xyoke, xcoil, xshoebase, xshoegap, xoffset, side, 'NWindingLayers', 3);

filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))


%% test any number of slots

slots = 2 * 21;
poles = 28;
ypole = 1;
ycoil = ypole / (slots / poles) * 0.8;
yshoegap = 0.2 * ycoil;
xyoke = 0.05;
xcoil = 0.3;
xshoebase = 0.05;
xshoegap = 0.025;
xoffset = 0;
side = 'l';

[FemmProblem, outernodes] = ...
    axialfluxstatorhalf2dfemmprob(slots, poles, ypole, ycoil, yshoegap, ...
                                  xyoke, xcoil, xshoebase, xshoegap, xoffset, side, ...
                                  'NWindingLayers', 2, ...
                                  'NSlots', 4);

filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))
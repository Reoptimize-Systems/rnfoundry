% Test_radialfluxstatorhalf2dfemmprob



Poles = 20;
slots = Poles*3;
thetapole = 2*pi/Poles;
thetacoil = (thetapole / (slots/Poles)) * 0.8;
thetashoegap = 0.2 * thetacoil;
ryoke = 0.02;
rcoil = 0.05;
rshoebase = 0.01;
rshoegap = 0.005;
roffset = 0.5;
side = 'i';

[FemmProblem, outernodes] = ...
    radialfluxstatorhalf2dfemmprob(slots, Poles, thetapole, thetacoil, thetashoegap, ...
                                   ryoke, rcoil, rshoebase, rshoegap, roffset, side, 'NWindingLayers', 1);

openprobleminfemm_mfemm( FemmProblem)



%% test any number of slots

slots = 2 * 21;
Poles = 28;
thetapole = 2*pi/Poles;
thetacoil = thetapole / (slots / Poles) * 0.8;
thetashoegap = 0.2 * thetacoil;
ryoke = 0.02;
rcoil = 0.05;
rshoebase = 0.01;
rshoegap = 0.005;
roffset = 0.5;
side = 'o';

[FemmProblem, outernodes] = ...
    axialfluxstatorhalf2dfemmprob(slots, Poles, thetapole, thetacoil, thetashoegap, ...
                                  ryoke, rcoil, rshoebase, rshoegap, roffset, side, ...
                                  'NWindingLayers', 2, ...
                                  'NSlots', 4);

filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))
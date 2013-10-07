% Test_axialfluxinnerstator2dfemmprob

slots = 2 * 21;
Poles = 28;
ypole = 1;
ycoil = ypole / (slots / Poles) * 0.8;
yshoegap = 0.2 * ycoil;
xyoke = 0.05;
xcoil = 0.3;
xshoebase = 0.05;
xshoegap = 0.025;


xmag = 0.1 * ypole;
xbackiron = 0.5 * xmag;
g = 0.05 * xmag;
yokecentresep = 2 * (xshoebase + xcoil + g + xbackiron + xmag) + xyoke;

2 * xmag + 2 * xbackiron + 2 * g

[FemmProblem, outernodes] = ...
    axialfluxinnerstator2dfemmprob(yokecentresep, slots, Poles, ypole, ...
        ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, 'NStators', 2);
    
    
filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))  

%% Any number of slots

slots = 2 * 21;
Poles = 28;
ypole = 1;
ycoil = ypole / (slots / Poles) * 0.8;
yshoegap = 0.2 * ycoil;
xyoke = 0.05;
xcoil = 0.3;
xshoebase = 0.05;
xshoegap = 0.025;


xmag = 0.1 * ypole;
xbackiron = 0.5 * xmag;
g = 0.05 * xmag;
yokecentresep = 2 * (xshoebase + xcoil + g + xbackiron + xmag) + xyoke;

2 * xmag + 2 * xbackiron + 2 * g

[FemmProblem, outernodes] = ...
    axialfluxinnerstator2dfemmprob(yokecentresep, slots, Poles, ypole, ...
        ycoil, yshoegap, xyoke, xcoil, xshoebase, xshoegap, 'NStators', 2, 'NSlots', 7);
    
    
filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))  


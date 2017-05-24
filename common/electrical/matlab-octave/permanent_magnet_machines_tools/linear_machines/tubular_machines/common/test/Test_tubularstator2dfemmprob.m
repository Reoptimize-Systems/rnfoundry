% Test_tubularstator2dfemmprob
%
% tubularstator2dfemmprob
%
%

slots = 6;
Poles = 2;
zpole = 1.0;
zcoil = 0.9 * 1/3;
zshoegap =  0.25;
ryoke = 0.2;
rcoil = 0.3;
rshoebase = 0.05;
rshoegap = rshoebase / 2;
roffset = 0.4;


[FemmProblem, info] = tubularstator2dfemmprob(slots, Poles, zpole, zcoil, zshoegap, ryoke, rcoil, rshoebase, rshoegap, roffset, 'l');

plotfemmproblem (FemmProblem)
% Test_slottedfemmprob_radial
%
%

design.StatorType = 'so';
design.Poles = 28;
design.Phases = 3;
design.Qs = design.Phases * 2 * 7;
design.yd = 1;
design.thetap = 2*pi/design.Poles;
design.thetam = design.thetap * 0.8;
design.thetac = (2*pi / design.Qs) * 0.8;
design.thetasg = design.thetac * 0.5;
design.tm = 0.01;
design.tbi = 0.01;
design.ty = 0.01;
design.tc = 0.03;
design.tsb = 0.01;
design.tsg = 0.005;
design.g = 3/1000;
design.ls = 0.3;

if strcmp(design.StatorType, 'si')
    design.Rmo = 0.5;
    design.Rmi = design.Rmi - design.tm;
    design.Rmm = mean([design.Rmi, design.Rmo]);
    design.Rci = design.Rmo + design.g + design.tsb;
    design.Rco = design.Rci + design.tc;
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbo = design.Rmi;
    design.Rbi = design.Rbo - design.tbi;
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Ryi = design.Rco;
    design.Ryo = design.Rco + design.ty;
    design.Rym = mean([design.Ryi, design.Ryo]);
elseif strcmp(design.StatorType, 'so')
    design.Rmi = 0.5;
    design.Rmo = design.Rmi + design.tm;
    design.Rmm = mean([design.Rmi, design.Rmo]);
    design.Rco = design.Rmi - design.g - design.tsb;
    design.Rci = design.Rco - design.tc;
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbi = design.Rmo;
    design.Rbo = design.Rbi + design.tbi;
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Ryo = design.Rci;
    design.Ryi = design.Ryo - design.ty;
    design.Rym = mean([design.Ryi, design.Ryo]);
end
%%
[FemmProblem, outermagsep, coillabellocs, yokenodeids] = ...
    slottedfemmprob_radial(design, ...
                           'StatorType', design.StatorType );

openprobleminfemm_mfemm(FemmProblem);

%%
[FemmProblem, coillabellocs] = slottedLfemmprob_radial(design, 'StatorType', design.StatorType);

openprobleminfemm_mfemm(FemmProblem);
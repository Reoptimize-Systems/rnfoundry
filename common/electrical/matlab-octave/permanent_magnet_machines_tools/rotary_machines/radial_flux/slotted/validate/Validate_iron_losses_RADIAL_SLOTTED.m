% Validate_iron_losses_RADIAL_SLOTTED
%
% losses will be validated by comparison with results taken from the paper:
%
% MI et al., "MODELING OF IRON LOSSES OF PM SYNCHRONOUS MOTORS", IEEE
% TRANSACTIONS ON INDUSTRY APPLICATIONS, VOL. 39, NO. 3, MAY/JUNE 2003
%
%

outputdir = '/home/rcrozier/Dropbox--Nova-Innovation/KTP/Outputs/iron-losses-validation-IIII'

%% set up design parameters
    
clear design simoptions
    
design.Poles = 4;
design.Qs = 36;
design.CoilLayers = 2;
design.Qc = design.Qs;
design.Phases = 3;
design.yd = 4;

design.Ryo = 95e-3;
design.g = 2e-3;
design.ty = 17.4e-3;
design.tm = 6.3e-3;
design.Rai = 58.5e-3;
design.Rmo = design.Rai - design.g;
design.tbi = 29.9523e-3;
design.tsb = 2.604e-3;
design.tsg = 1.7364e-3;
design.tc = 19.1e-3 - design.tsb;
design.Rtsb = design.Rai + design.tsb;
design.Ryi = design.Rtsb + design.tc(1);
design.thetam = (tau / design.Poles) * 0.667;
design.tc(2) = design.tc(1) * 0.133333;

chordanglefcn = @(r,l) 2 .* asin(l / (2.*r));

design.thetacg = feval(chordanglefcn, design.Rtsb, 5.178e-3);
design.thetacy = feval(chordanglefcn, design.Ryi, 7.2156e-3);
design.thetasg = feval(chordanglefcn, design.Ryi, 3e-3);
design.ls = 76e-3; % paper value: 88.9e-3;

design.ArmatureType = 'external';
design.MagnetPolarisation = 'radial';

design = completedesign_RADIAL_SLOTTED(design);

design.MagFEASimMaterials.Magnet = matstr2matstruct_mfemm ( 'NdFeB 40 MGOe' );
design.MagFEASimMaterials.Magnet.H_c = 815000;

design.MagFEASimMaterials.FieldBackIron = 'M-19 Steel';
design.MagFEASimMaterials.ArmatureYoke = design.MagFEASimMaterials.FieldBackIron;
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';
design.MagFEASimMaterials.AirGap = 'Air';
    
%% get the machine data

design.CoilFillFactor = 0.6;
design.CoilTurns = 100;
design.RlVRp = 100;

simoptions = struct();
simoptions.GetVariableGapForce = false;
simoptions.NMagFEAPositions = 30;
simoptions.MagFEASim.BackIronRegionMeshSize = ...
    choosemesharea_mfemm(design.ty, design.Rai*design.thetap, 1/30);

simoptions.reltol = 1e-4;
%simoptions.PhaseCurrentTols = repmat(0.001, 1, design.Phases);
%simoptions.maxstep = (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1)) / 1000;
simoptions.PoleCount = 1000;

design.CoreLoss = struct ();

[design, simoptions] = simfun_RADIAL_SLOTTED (design, simoptions);
simoptions.ODESim.PreProcFcn = [];

save (fullfile(outputdir, 'validate_iron_losses_design_and_simoptions.mat'), 'design', 'simoptions');

% MiIronLosses = [];

rpm = [ 300, 600, 900, 1200, 1500, 1800 ];
vel = rpm2vel (rpm, design.Rmo + design.g/2);
freq = rpm .* 2 ./ 60;
experim = [ 3.1; 11.9; 21.1; 35.4; 55.0; 72.5 ];
MiIronLosses = [ 6.7; 15.6; 26.7; 40.1; 55.8; 73.7 ];

colheadings = {  '', 'm/s' ,  'Hz' , '' , 'Mi Model' , 'Mi Coeff', 'M-19 29', 'M-19 26' , 'M-19 24' , 'M-36 29' , 'M-36 26' , 'M-45 29'};
wid = 14;
fms = {'.4g'};

colsep = ' & ';
rowending = ' \\';

fileID = 1;

%% Using Mi Coefficients

design.CoreLoss.kh =  44;
design.CoreLoss.kc =  0.07;
design.CoreLoss.ke =  4.0267e-11;
design.CoreLoss.beta =  2;
ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
MiIronLosses = [ MiIronLosses, ptables.PowerLossIronMean ];


% design.CoreLoss.kh =  44 * 2.2;
% design.CoreLoss.kc =  0.07 * 2.2;
% design.CoreLoss.ke =  4.0267e-11 * 2.2;
% design.CoreLoss.beta =  2;
% ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
% MiIronLosses = [ MiIronLosses, ptables.PowerLossIronMean ];

%% AK Steel material data
simpleironlosses = MiIronLosses;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
                                      
export_fig (fullfile(outputdir, 'M-19-29-freqrange-full-Brange-full.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
simpleironlosses = [ simpleironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '26', ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
                                      
export_fig (fullfile(outputdir, 'M-19-26-freqrange-full-Brange-full.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
simpleironlosses = [ simpleironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '24', ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-19-24-freqrange-full-Brange-full.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
simpleironlosses = [ simpleironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '29', ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-36-29-freqrange-full-Brange-full.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
simpleironlosses = [ simpleironlosses, ptables.PowerLossIronMean ];

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-36-26-freqrange-full-Brange-full.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
simpleironlosses = [ simpleironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-45', '29', ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-45-29-freqrange-full-Brange-full.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
simpleironlosses = [ simpleironlosses, ptables.PowerLossIronMean ];

modelerror = 100 * bsxfun(@rdivide, bsxfun(@minus, simpleironlosses, experim), experim);

displaytable([rpm', vel', freq', experim, simpleironlosses(:,1:2), simpleironlosses(:,3:end)],colheadings,wid,fms,{},fileID,colsep,rowending);
fprintf (1, '\\midrule\n');
displaytable([rpm', vel', freq', zeros(size(modelerror,1),1), modelerror(:,1:2), modelerror(:,3:end)],{},wid,fms,{},fileID,colsep,rowending);

save (fullfile(outputdir, 'validate_iron_lossses_intermediate.mat'));

%% With added zeros

zerosironlosses = MiIronLosses;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'InterpolateMissing', false, 'AddZeros', true, 'Plot', '3DLog');

export_fig (fullfile(outputdir, 'M-19-29-freqrange-full-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zerosironlosses = [ zerosironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '26', 'AddZeros', true, 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-19-26-freqrange-full-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zerosironlosses = [ zerosironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '24', 'AddZeros', true, 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-19-24-freqrange-full-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zerosironlosses = [ zerosironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '29', 'AddZeros', true, 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-36-29-freqrange-full-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zerosironlosses = [ zerosironlosses, ptables.PowerLossIronMean ];

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'AddZeros', true, 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-36-26-freqrange-full-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zerosironlosses = [ zerosironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-45', '29', 'AddZeros', true, 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-45-29-freqrange-full-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zerosironlosses = [ zerosironlosses, ptables.PowerLossIronMean ];


save (fullfile(outputdir, 'validate_iron_lossses_intermediate.mat'));

modelerror = 100 * bsxfun(@rdivide, bsxfun(@minus, zerosironlosses, experim), experim);

displaytable([rpm', vel', freq', experim, zerosironlosses(:,1:2), zerosironlosses(:,3:end)],colheadings,wid,fms,{},fileID,colsep,rowending);
fprintf (1, '\\midrule\n');
displaytable([rpm', vel', freq', zeros(size(modelerror,1),1), modelerror(:,1:2), modelerror(:,3:end)],{},wid,fms,{},fileID,colsep,rowending);
                                                     

                                                     
%% With added zeros and limited range of frequencies

zeroslimitironlosses = MiIronLosses;
limitfreqs = 600;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');

export_fig (fullfile(outputdir, 'M-19-29-freqrange-0-600-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '26', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-19-26-freqrange-0-600-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '24', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-19-24-freqrange-0-600-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-36-29-freqrange-0-600-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-36-26-freqrange-0-600-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-45', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, 'M-45-29-freqrange-0-600-Brange-full-zeros.pdf'), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

save (fullfile(outputdir, 'validate_iron_lossses_intermediate.mat'));

modelerror = 100 * bsxfun(@rdivide, bsxfun(@minus, zeroslimitironlosses, experim), experim);

fprintf (1, '%d Positons, Freqs limited to %g Hz\n\n', simoptions.NMagFEAPositions, limitfreqs);
displaytable([rpm', vel', freq', experim, zeroslimitironlosses(:,1:2), zeroslimitironlosses(:,3:end)],colheadings,wid,fms,{},fileID,colsep,rowending);
fprintf (1, '\\midrule\n');
displaytable([rpm', vel', freq', zeros(size(modelerror,1),1), modelerror(:,1:2), modelerror(:,3:end)],{},wid,fms,{},fileID,colsep,rowending);

%%
zeroslimitironlosses = MiIronLosses;
limitfreqs = 400;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-29-freqrange-0-', int2str(limitfreqs), '-Brange-full-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '26', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-26-freqrange-0-', int2str(limitfreqs), '-Brange-full-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '24', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-24-freqrange-0-', int2str(limitfreqs), '-Brange-full-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-36-29-freqrange-0-', int2str(limitfreqs), '-Brange-full-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-36-26-freqrange-0-', int2str(limitfreqs), '-Brange-full-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-45', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-45-29-freqrange-0-', int2str(limitfreqs), '-Brange-full-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

save (fullfile(outputdir, 'validate_iron_lossses_intermediate.mat'));

modelerror = 100 * bsxfun(@rdivide, bsxfun(@minus, zeroslimitironlosses, experim), experim);

fprintf (1, '%d Positons, Freqs limited to %g Hz\n\n', simoptions.NMagFEAPositions, limitfreqs);
displaytable([rpm', vel', freq', experim, zeroslimitironlosses(:,1:2), zeroslimitironlosses(:,3:end)],colheadings,wid,fms,{},fileID,colsep,rowending);
fprintf (1, '\\midrule\n');
displaytable([rpm', vel', freq', zeros(size(modelerror,1),1), modelerror(:,1:2), modelerror(:,3:end)],{},wid,fms,{},fileID,colsep,rowending);

eval (['zeros_limit_freq_ironlosses_', int2str(limitfreqs), ' = zeroslimitironlosses;']);

%%
zeroslimitironlosses = MiIronLosses;
limitB = 1.2;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'AddZeros', true, ...
                                          'LimitInduction', [0, limitB], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-29-freqrange-full-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '26', 'AddZeros', true, ...
                                          'LimitInduction', [0, limitB], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-26-freqrange-full-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '24', 'AddZeros', true, ...
                                          'LimitInduction', [0, limitB], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-24-freqrange-full-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '29', 'AddZeros', true, ...
                                          'LimitInduction', [0, limitB], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-36-29-freqrange-full-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'AddZeros', true, ...
                                          'LimitInduction', [0, limitB], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-36-26-freqrange-full-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-45', '29', 'AddZeros', true, ...
                                          'LimitInduction', [0, limitB], 'InterpolateMissing', false, 'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-45-29-freqrange-full-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

save (fullfile(outputdir, 'validate_iron_lossses_intermediate.mat'));

modelerror = 100 * bsxfun(@rdivide, bsxfun(@minus, zeroslimitironlosses, experim), experim);

fprintf (1, '%d Positons, Induction limited to %g T\n\n', simoptions.NMagFEAPositions, limitB);
displaytable([rpm', vel', freq', experim, zeroslimitironlosses(:,1:2), zeroslimitironlosses(:,3:end)],colheadings,wid,fms,{},fileID,colsep,rowending);
fprintf (1, '\\midrule\n');
displaytable([rpm', vel', freq', zeros(size(modelerror,1),1), modelerror(:,1:2), modelerror(:,3:end)],{},wid,fms,{},fileID,colsep,rowending);

eval (['zeros_limitinduc_ironlosses_', int2str(1000*limitB), ' = zeroslimitironlosses;']);


%%
zeroslimitironlosses = MiIronLosses;
limitB = 1.2;
limitfreqs = 600;

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], ...
                                          'LimitInduction', [0, limitB], ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-29-freqrange-0-',int2str(limitfreqs),'-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '26', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], ...
                                          'LimitInduction', [0, limitB], ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-26-freqrange-0-',int2str(limitfreqs),'-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-19', '24', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], ...
                                          'LimitInduction', [0, limitB], ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-19-24-freqrange-0-',int2str(limitfreqs),'-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], ...
                                          'LimitInduction', [0, limitB], ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-36-29-freqrange-0-',int2str(limitfreqs),'-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];

[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], ...
                                          'LimitInduction', [0, limitB], ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-36-26-freqrange-0-',int2str(limitfreqs),'-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];


[design.CoreLoss.kh, ...
 design.CoreLoss.kc, ...
 design.CoreLoss.ke, ...
 design.CoreLoss.beta ] = corelosscoeffs ('M-45', '29', 'AddZeros', true, ...
                                          'LimitFreqs', [0, limitfreqs], ...
                                          'LimitInduction', [0, limitB], ...
                                          'InterpolateMissing', false, ...
                                          'Plot', '3DLog');
 
export_fig (fullfile(outputdir, ['M-45-29-freqrange-0-',int2str(limitfreqs),'-Brange-0-',int2str(1000*limitB),'-zeros.pdf']), '-pdf');
close (gcf)

ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, design.RlVRp);
zeroslimitironlosses = [ zeroslimitironlosses, ptables.PowerLossIronMean ];



modelerror = 100 * bsxfun(@rdivide, bsxfun(@minus, zeroslimitironlosses, experim), experim);

fprintf (1, '%d Positons, Induction limited to %g B and %g Hz\n\n', simoptions.NMagFEAPositions, limitB, limitfreqs);
displaytable([rpm', vel', freq', experim, zeroslimitironlosses(:,1:2), zeroslimitironlosses(:,3:end)],colheadings,wid,fms,{},fileID,colsep,rowending);
fprintf (1, '\\midrule\n');
displaytable([rpm', vel', freq', zeros(size(modelerror,1),1), modelerror(:,1:2), modelerror(:,3:end)],{},wid,fms,{},fileID,colsep,rowending);

save (fullfile(outputdir, 'validate_iron_lossses_intermediate.mat'));

eval (['zeros_limitinduc_freq_ironlosses_', int2str(1000*limitB), '_', int2str(1000*limitfreqs), ' = zeroslimitironlosses;']);

%% save the data

save (fullfile(outputdir, 'validate_iron_losses_RADIAL_SLOTTED.mat'));



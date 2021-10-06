% Test_slottedfemmprob_radial
%
%

design = test_design_RADIAL_SLOTTED ();

design = completedesign_RADIAL_SLOTTED (design);

    
%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'DrawCoilInsulation', true, ...
                             'NPolePairs', 2);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%% Full number of machine poles (but not 'Full' DrawingType)

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'DrawCoilInsulation', true, ...
                             'NPolePairs', 34);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

% %%
% 
% [FemmProblem, coillabellocs] = ...
%     slottedLfemmprob_radial (design, ...
%                              'ArmatureType', design.ArmatureType );
% 
% plotfemmproblem (FemmProblem);
% % openprobleminfemm_mfemm (FemmProblem);
% 
% %%
% [FemmProblem, coillabellocs] = slottedLfemmprob_radial (design, 'ArmatureType', design.ArmatureType);
% 
% openprobleminfemm_mfemm (FemmProblem);

%%

design.ShoeCurveControlFrac = 0.9;
design.CoilInsulationThickness = 2/10000;

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'DrawOuterRegions', false );
                             
openprobleminfemm_mfemm (FemmProblem);


%%

design = test_design_RADIAL_SLOTTED ();
design = completedesign_RADIAL_SLOTTED (design);

design.qc = fr (1,4);
design.yd = 1;
design.NBasicWindings = 10;
design.CoilLayers = 1;
design = rmfield (design, {'Qs', 'Qc', 'Poles', 'g', 'thetap', 'thetam', 'Ryi', 'Rai', 'Rmo', 'Rmi'});
design = completedesign_RADIAL_SLOTTED (design);

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'DrawCoilInsulation', false, ...
                             'DrawOuterRegions', true );
                             
openprobleminfemm_mfemm (FemmProblem);


%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', 10, ...
                             'DrawingType', 'BoundarySwap');

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', 1, ...
                             'DrawingType', 'BoundarySwap', ...
                             'BoundaryShift', 4);

% plotfemmproblem (FemmProblem);
openprobleminfemm_mfemm (FemmProblem);

%% No tsg, tsb, thetasg etc.

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design = rmiffield (design, 'tsb');
design = rmiffield (design, 'tsg');
design = rmiffield (design, 'thetasg');

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', 2);

plotfemmproblem (FemmProblem);

%% No base curve at all

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.tc(2) = 0;

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', 2);

plotfemmproblem (FemmProblem);

%% No base curve at all, full circumference

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.tc(2) = 0;

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', design.Poles/2);

plotfemmproblem (FemmProblem);

%% Any number of external regions

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';


[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', 2, ...
                             'StatorOuterRegionSize', [design.tm, 3*design.tm, design.tm, 2*design.tm] );

plotfemmproblem (FemmProblem);

%% Any number of external regions with specified materials in the external 
%% regions

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.tc(2) = 0;

[FemmProblem, coillabellocs] = ...
    slotlessfemmprob_radial ( design, ...
                             'NPolePairs', 2, ...
                             'StatorOuterRegionSize', [design.tm, 3*design.tm, design.tm, 2*design.tm], ...
                             'StatorOuterRegionMaterials', { design.MagFEASimMaterials.ArmatureCoil, ...  
                                                             design.MagFEASimMaterials.ArmatureYoke, ...
                                                             design.MagFEASimMaterials.CoilInsulation, ...
                                                             design.MagFEASimMaterials.FieldBackIron } );

plotfemmproblem (FemmProblem);


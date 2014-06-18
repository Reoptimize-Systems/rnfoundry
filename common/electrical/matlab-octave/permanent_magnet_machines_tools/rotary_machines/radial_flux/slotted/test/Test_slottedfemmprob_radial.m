% Test_slottedfemmprob_radial
%
%

design = test_design_RADIAL_SLOTTED ();

design = completedesign_RADIAL_SLOTTED(design);

    
%%
[FemmProblem, coillabellocs, yokenodeids] = ...
    slottedfemmprob_radial(design, ...
                           'ArmatureType', design.ArmatureType );

openprobleminfemm_mfemm(FemmProblem);

[FemmProblem, coillabellocs, yokenodeids] = ...
    slottedLfemmprob_radial(design, ...
                           'ArmatureType', design.ArmatureType );

openprobleminfemm_mfemm(FemmProblem);

%%
[FemmProblem, coillabellocs] = slottedLfemmprob_radial(design, 'ArmatureType', design.ArmatureType);

openprobleminfemm_mfemm(FemmProblem);

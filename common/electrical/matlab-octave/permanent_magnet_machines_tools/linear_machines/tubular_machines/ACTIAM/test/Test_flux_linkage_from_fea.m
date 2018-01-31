
% Test_feaprescribedmotodeforcefcn_linear.m

clear design simoptions

% design = test_design_TM_SLOTLESS ('dims');

load ('/home/rcrozier/Sync/work/enercro/Projects/WES-Wavedrive/Outputs/design_004_converted_to_rnfoundry.mat')

%%
design.RlVRp = 8;

% [design.CoreLoss.kh, ...
%  design.CoreLoss.kc, ...
%  design.CoreLoss.ke, ...
%  design.CoreLoss.beta ] = corelosscoeffs ('M-19', '29', 'InterpolateMissing', false);


simoptions = struct ();
simoptions.GetVariableGapForce = false;
simoptions.MagFEASimType = 'multiple';
simoptions.UseParFor = true;
simoptions.MagFEASim.UseParFor = false;

design = completedesign_TM_SLOTLESS (design, simoptions);

design.zs = design.zp / 3;

% setup simulation options
%simoptions = simsetup_ROTARY(design, 'simfun_RADIAL_SLOTTED', 'finfun_RADIAL_SLOTTED', ...
%                                'Rpm', 30, ...
%                                'TSpan', [0,1]);

simoptions = simsetup_linear (design, 'simfun_TM_SLOTLESS', 'feaprescribedmotfinfun_TM_SLOTLESS', ...
                                'ForceFcn', 'forcefcn_linear_pscbmot', ...
                                'EvalFcn', 'feaprescribedmotodeforcefcn_linear', ...
                                'simoptions', simoptions, ...
                                'PoleCount', 30, ...
                                ... 'RampPoles', 5, ...
                                'Velocity', 1.0 );
                            
simoptions.ODESim.MaxStep = simoptions.ODESim.MaxStep / 10;

% simoptions.ODESim.OutputFcn = 'feaprescribedmotodeoutputfcn_ROTARY';

design.FEAFluxLinkageFCN = @feaode_TM_SLOTLESS;

simoptions.reltol = 1e-6;
simoptions.PhaseCurrentTols = repmat (0.001, 1, design.Phases);
% simoptions.maxstep = (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1)) / 10000;

simoptions.Evaluation = designandevaloptions_TM_SLOTLESS ();

[design, simoptions] = feval (simoptions.ODESim.PreProcFcn, design, simoptions);
simoptions.ODESim.PreProcFcn = [];

[design, simoptions] = feval (simoptions.ODESim.PostPreProcFcn, design, simoptions);
simoptions.ODESim.PostPreProcFcn = [];


simoptions.UseParFor = false;
simoptions.MagFEASim.UseParFor = false;

%%

npoles = 3;
feapos = linspace (-npoles, npoles, 2*npoles * 20);

[ RawCoggingForce, ...
  BxCoreLossData, ...
  ByCoreLossData, ...
  ArmatureToothFluxDensity, ...
  FemmDirectFluxLinkage, ...
  AslotPos, ...
  slotIntA, ...
  BslotPos, ...
  slotIntB, ...
  design ] = feasim_TM_SLOTLESS (design, simoptions, design.zp * feapos(1), ...
                                    'IsInitialisation', true, ...
                                    'NPolePairs', 1);

simoptions.MagFEASim = setfieldifabsent ( simoptions.MagFEASim, 'UseParFor', false);

if simoptions.MagFEASim.UseParFor

    parfor posind = 2:numel (feapos)

        [ ~, ...
          ~, ...
          ~, ...
          ~, ...
          FemmDirectFluxLinkage(posind,:)] = feasim_TM_SLOTLESS (design, simoptions, design.zp * feapos(posind), ...
                                                'IsInitialisation', false, ...
                                                'NPolePairs', 1);

    end

else

    for posind = 2:numel (feapos)

        [ ~, ...
          ~, ...
          ~, ...
          ~, ...
          FemmDirectFluxLinkage(posind,:)] = feasim_TM_SLOTLESS (design, simoptions, design.zp * feapos(posind), ...
                                                'IsInitialisation', false, ...
                                                'NPolePairs', 1);

    end

end

%%

xE = design.zp/2 + 2*design.zp/3;

plotyy (feapos, FemmDirectFluxLinkage(:,1), ...
    feapos, periodicslmeval (feapos + xE / design.zp, design.slm_fluxlinkage));

% legend ('1', '2', '3', '4');

%%

for posind = 1:numel (feapos)
    [~, ~, ~, dlambdaVdt(posind,:), ~, ~] = machineodesim_AM(design, simoptions, design.zp*feapos(posind) + xE, 0, 1.0, 0, [0;0;0]);
end

plotyy (feapos, FemmDirectFluxLinkage(:,1), ...
    feapos, dlambdaVdt(:,1));

%%

  [ ~, ...
          ~, ...
          ~, ...
          ~, ...
          FemmDirectFluxLinkage(posind,:)] = feasim_TM_SLOTLESS (design, simoptions, design.zp * feapos(22), ...
                                                'IsInitialisation', false, ...
                                                'NPolePairs', 1);
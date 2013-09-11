%% Test_getfeadata_PMSM
%
% A script for testing the function getfeadata_PMSM
%
%
% First we need some plausible machine variables for testing, we will use
% the AWS variables for this

clear design simoptions

design.phases = 3;
design.Wp = 0.12;
design.Wm = 0.8*design.Wp;
design.lm = 0.015;

design.Hc = 979000;

design.fillfactor = 0.585;
design.g = 0.003; 
design.ls = 1; 
design.Dc = 0.005;

design.hba = design.Wp / 4;
design.hbf = design.Wp / 4; 
design.Ws = design.Wp / design.phases; 
design.Wt = design.Ws / 2;
design.ht = 5*design.Wt;
design.hm = 0.2 * design.Wm;

% Number of turns
design.CoilTurns = 200;

design.FEMMTol = 1e-5;

% Get the dimensionless ratios from the parameters
design = dimensions2ratios_PMSM(design);

% set design mode
design.mode = [1, 1, 0, 1];

Jcoil = 0;
xRVTaup = linspace(0,2,50);

if ~isfemmopen
    openfemm;
end

design = setfieldifabsent(design, 'CoilLayers', 2);
    
if design.CoilLayers == 1
    design.mode(4) = 0;
elseif design.CoilLayers == 2
    design.mode(4) = 1;
else
    error('Unrecognised number of coil layers');
end

% set up variables for calculation of the iron losses, if they have not
% been supplied
if ~isfield(design, 'CoreLoss')
    % Coreloss(1) will be the armature core data
    [design.CoreLoss.fq, ...
     design.CoreLoss.Bq, ...
     design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
    % CoreLoss(2) will be the translator core data, by default the same
    % as the armature
    design.CoreLoss(2) = design.CoreLoss(1);
end

[design, simoptions] = simfun_linear(design, struct());

% Determine the area of the coil, and the find either the number of
% turns, coil wire diameter if required
design = checkcoilprops_AM(design);

design.DcVWc = design.Dc ./ design.Wc;

% set up the positions at which we will gather data. These are shifted
% by 0.5 to get the flux linkage from it's maximum
    
[output, design] = getfeadata_PMSM(design, Jcoil, xRVTaup);

[yokehistloss, yokeeddyloss, yokeexcessloss] = ...
        softferrolossrectregionvarxpartcalc( design.CoreLoss(1).Bx, ...
                                             design.CoreLoss(1).By, ...
                                             design.CoreLoss(1).Bz, ...
                                             design.CoreLoss(1).Hx, ...
                                             design.CoreLoss(1).Hy, ...
                                             design.CoreLoss(1).Hz, ...
                                             design.CoreLoss(1).kc, ...
                                             design.CoreLoss(1).ke, ...
                                             design.CoreLoss(1).beta, ...
                                             design.CoreLoss(1).xstep, ...
                                             design.CoreLoss(1).dx, ...
                                             design.CoreLoss(1).dy, ...
                                             design.CoreLoss(1).dz );
                                         
[teethhistloss, teetheddyloss, teethexcessloss] = ...
        softferrolossrectregionvarxpartcalc( design.CoreLoss(2).Bx, ...
                                             design.CoreLoss(2).By, ...
                                             design.CoreLoss(2).Bz, ...
                                             design.CoreLoss(2).Hx, ...
                                             design.CoreLoss(2).Hy, ...
                                             design.CoreLoss(2).Hz, ...
                                             design.CoreLoss(2).kc, ...
                                             design.CoreLoss(2).ke, ...
                                             design.CoreLoss(2).beta, ...
                                             design.CoreLoss(2).xstep, ...
                                             design.CoreLoss(2).dx, ...
                                             design.CoreLoss(2).dy, ...
                                             design.CoreLoss(2).dz ); 
                                         
plot(xRVTaup, yokehistloss + teethhistloss);

slmengine(xRVTaup(xRVTaup <= 1), yokehistloss(xRVTaup <= 1) + teethhistloss(xRVTaup <= 1), 'plot', 'on', 'knots', 26, 'endcon', 'periodic')


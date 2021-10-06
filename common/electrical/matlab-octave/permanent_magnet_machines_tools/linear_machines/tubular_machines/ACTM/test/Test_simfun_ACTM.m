% Test_simfun_ACTM

clear design simoptions

design.Phases = 3;         % Number of Phases in machine
design.Rm = 0.1;
design.g = 5/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.5;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.025;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.CoilFillFactor = 0.65;
%design.Dc = 1/1000;  % 1 mm diameter wire 
design.CoilTurns = 500;
design.mode = 2; 
design.LgVLc = 0;
design.Poles = [10 30];
% design.FieldDirection = 1;
% design.PowerPoles = Poles(1);

design = ratios2dimensions_ACTM(design);

[design, simoptions] = simfun_ACTM(design, struct());


plotfemmproblem(design.FemmProblem);

hold on
contourf(design.X, design.Y .* design.Wp, design.A);
hold off

figure; scatter3(design.X(:), design.Y(:), design.Bx(:));

figure; scatter3(design.X(:), design.Y(:), design.By(:));

figure; scatter3(design.X(:), design.Y(:), design.A(:));


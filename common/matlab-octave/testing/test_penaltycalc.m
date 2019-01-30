%% Upper
vars.bananas = 10
penaltyspec.bananas.Limit = 8
penaltyspec.bananas.PenaltyFactors = 5
penaltyspec.bananas.Type = 'upper';

penalties = opt.penaltycalc (vars, penaltyspec)

%% Lower
vars.bananas = 5
penaltyspec.bananas.Limit = 8
penaltyspec.bananas.PenaltyFactors = 5
penaltyspec.bananas.Type = 'lower';

penalties = opt.penaltycalc (vars, penaltyspec)



%% Target

vars.bananas = 6
penaltyspec.bananas.Limit = 8
penaltyspec.bananas.PenaltyFactors = 5
penaltyspec.bananas.Type = 'target';

penalties = opt.penaltycalc (vars, penaltyspec)


vars.bananas = 10
penaltyspec.bananas.Limit = 8
penaltyspec.bananas.PenaltyFactors = 5
penaltyspec.bananas.Type = 'target';

penalties = opt.penaltycalc (vars, penaltyspec)
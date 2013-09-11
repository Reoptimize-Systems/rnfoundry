function tf = isslm(teststruct)
% isslm: test for structure fields indicating a structure is an slm object

% Created by Richard Crozier 2012

    tf = every(isfield(teststruct, {'form', 'degree', 'knots', 'coef', 'stats'}));

end
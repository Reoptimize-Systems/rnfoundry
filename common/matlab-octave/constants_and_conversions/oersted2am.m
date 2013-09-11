function ham = oersted2am(hoe)
% converts magnetizing field strength in Oersted c.g.s units to
% Amperes/metre in S.I. units

% Created by Richard Crozier 

    ham = hoe * 1000 / (4 * pi);

end
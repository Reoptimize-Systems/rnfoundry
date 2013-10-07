% Test_febeamdef_FM

clear 
% cross-sectional dimensions and length of the beams in metres
b = 0.05; 
h = 0.10; 

% Table20r4K([h,b])

if b >= h
    a1 = b/2;
    a2 = h/2;
else    
    a1 = h/2;
    a2 = b/2;
end

% at ends of field magnet supports
BeamInfo.SupportBeams = struct('A', b * h, ...
                               'I2', b * h^3 / 12, ...
                               'I3', b^3 * h / 12, ...
                               'I1', b^3 * h / 3, ...
                               'J',  (a1)*(a2)^3 * ((16/3) - 3.36*(a2/a1)*(1 - (a2^4 / (12*a1^4)))), ...
                               'rho', 7500, ...
                               'AngleFromHorizontal', pi/2);  
                           
% % at ends of field magnet supports
% BeamInfo.SupportBeams = struct('A', b * h, ...
%                                'I2', b * h^3 / 12, ...
%                                'I3', b^3 * h / 12, ...
%                                'I1', b^3 * h / 3);                             

fprintf(1, '\nI1: %e, J: %e\n', BeamInfo.SupportBeams.I1, BeamInfo.SupportBeams.J)

%Young's modulus of the beam material in Pa 
BeamInfo.SupportBeams.E = 207e9;

% Poisson ratio
BeamInfo.SupportBeams.nu = 0.31;

BeamInfo.GuideRails = BeamInfo.SupportBeams;
BeamInfo.GuideRails.length = 12;
BeamInfo.GuideRails.Sections = 5;

BeamInfo.GuideBearings = BeamInfo.SupportBeams;
BeamInfo.GuideBearings.NoPerGuide = 2;
BeamInfo.GuideBearings.E = BeamInfo.GuideBearings.E * 100;

% holds the field sides apart
BeamInfo.OuterWebs = BeamInfo.SupportBeams;
BeamInfo.OuterWebs.NoPerSide = 2;
BeamInfo.OuterWebs.Sections = 3;

% field magnet supports
BeamInfo.OuterPoleSupports = BeamInfo.SupportBeams;
BeamInfo.OuterPoleSupports.Sections = 10;
BeamInfo.OuterPoleSupports.NoPerSide = 8;

% x = design.ls .* options.alphab ./ 2;
% 
% y = (design.hbf + design.hm + design.g + desing.ht + design.hba);
% The height of the frame above zero is half the total translator height
% z = ((design.Poles(2) .* design.Wp) - (design.Wp/2)) / 2;

BeamInfo.width = 2;
BeamInfo.depth = 1.0;
BeamInfo.height = 6.0;

MForce = 1000;

tic
[airgapdisps, BeamInfo2, fens, gcells] = febeamdef_FM(BeamInfo, MForce, true);
toc
%%

[airgapdisps, BeamInfo2, fens, gcells] = febeamdef_FM(BeamInfo2, MForce, false, fens, gcells);


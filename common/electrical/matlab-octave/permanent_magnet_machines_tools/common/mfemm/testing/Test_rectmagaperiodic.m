% Test_planarrectmagaperiodic

FemmProblem = newproblem_mfemm('planar');

%   Materials
Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

FemmProblem.Materials = Matlib([1, 47, 2]);

ypole = 1; 
ymag = 0.8;
xmag = 0.2;
xoffset = -4;

%% cases

test = 4;

switch test

    case 1
        % no displacement
        pos = 0;

    case 2
        % B top on top
        pos = (ypole - ymag)/2;

    case 3
        % A bot on bottom
        pos = -(ypole - ymag)/2;

    case 4
        % halway B
        pos = (ymag/2 + (ypole - ymag)/2);

    case 5
        % halway A
        pos = -(ymag/2 + (ypole - ymag)/2);

    case 6
        % B bot on bot 
        pos = (ymag + (ypole - ymag)/2);

    case 7
        % A top on top
        pos = -(ymag + (ypole - ymag)/2);
        
    case 8
        % one period translation
        % A bot B top
        pos = 2*(ymag + (ypole - ymag));       
        
    case 9
        % -ve one period translation
        % A bot B top
        pos = -2*(ymag + (ypole - ymag));   
        
    case 10
        % half period translation
        % A top B bot
        pos = (ymag + (ypole - ymag));     
        
    case 11
        % -ve half period translation
        % A top B bot
        pos = -(ymag + (ypole - ymag));           
        
end

%%

[FemmProblem, nodes, links] = planarrectmagaperiodic(FemmProblem, ypole, ...
                                                     ymag, xmag, xoffset, ...
                                                     pos, 'MagMaterial', 3, ...
                                                     'MagDirections', [90, 270]);

% plotnodelinks(nodes, links)

filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))


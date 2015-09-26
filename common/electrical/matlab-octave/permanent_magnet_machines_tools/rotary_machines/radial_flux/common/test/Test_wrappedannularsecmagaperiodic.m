% Test_wrappedannularsecmagaperiodic


thetapole = 2*pi/20; 
thetamag = thetapole * 0.8;
rmag = 0.02;
roffset = 0.5;
wrapperthickness = [0.1, 0.2; 0.1, 0.4; 0.2, 0.5];


%% cases

test = 4;

switch test

    case 1
        % no displacement
        pos = 0;

    case 2
        % B top on top
        pos = (thetapole - thetamag)/2;

    case 3
        % A bot on bottom
        pos = -(thetapole - thetamag)/2;

    case 4
        % halway B
        pos = (thetamag/2 + (thetapole - thetamag)/2);

    case 5
        % halway A
        pos = -(thetamag/2 + (thetapole - thetamag)/2);

    case 6
        % B bot on bot 
        pos = (thetamag + (thetapole - thetamag)/2);

    case 7
        % A top on top
        pos = -(thetamag + (thetapole - thetamag)/2);
        
    case 8
        % one period translation
        % A bot B top
        pos = 2*(thetamag + (thetapole - thetamag));       
        
    case 9
        % -ve one period translation
        % A bot B top
        pos = -2*(thetamag + (thetapole - thetamag));   
        
    case 10
        % half period translation
        % A top B bot
        pos = (thetamag + (thetapole - thetamag));     
        
    case 11
        % -ve half period translation
        % A top B bot
        pos = -(thetamag + (thetapole - thetamag));           
        
end

%%

FemmProblem = newproblem_mfemm('planar');

%   Materials
Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

FemmProblem.Materials = Matlib([1, 47, 2]);


[FemmProblem, wrapperthickness, info] = wrappedannularsecmagaperiodic(FemmProblem, thetapole, ...
                                                          thetamag, rmag, roffset, ...
                                                          pos, wrapperthickness, ...
                                                          'MagnetMaterial', 3);

% plotnodelinks(nodes, links)

plotfemmproblem (FemmProblem);

%%

FemmProblem = newproblem_mfemm('planar');

%   Materials
Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

FemmProblem.Materials = Matlib([1, 47, 2]);


[FemmProblem, wrapperthickness, info] = wrappedannularsecmagaperiodic(FemmProblem, thetapole, ...
                                                          thetamag, rmag, roffset, ...
                                                          pos, wrapperthickness, ...
                                                          'MagnetMaterial', 3, ...
                                                          'NPolePairs', 2);

% plotnodelinks(nodes, links)

plotfemmproblem (FemmProblem);

%% Test_mbdyn_MBCNodal


% start mbdyn 
delete ('output.*');
[status, cmdout] = system (sprintf ('mbdyn -f %s -o output > output.txt 2>&1 &', fullfile (getmfilepath ('Test_mbdyn_MBCNodal'), 'socket.mbd')))

% wait a few seconds for mbdyn to initialise
pause (3);

%%

mb = mbdyn.MBCNodal ();

mb.Initialize ( 'local', '/tmp/mbdyn.sock', ...
                'NNodes', 2, ...
                'UseLabels', false);

nnodes = mb.GetNodes ()



maxiter = 3; % to match example program 

forces = [ 0., 0., 1.;
           0., 0., 1. ].';

status = 0;
while status == 0
    
    for iter = (maxiter-1):-1:1

        disp(iter)
        status = mb.GetMotion ()

        pos1 = mb.X(1)

        % set the forces
        mb.F (forces)

        mb.applyForcesAndMoments (false)

    end
    
    status = mb.GetMotion ()

    pos1 = mb.X(1)

    % set the forces
    mb.F (forces)

    mb.applyForcesAndMoments (true)

end

% 
% mb.GetMotion ();
% 
% pos1 = mb.X(1)
% 
% mb.F (forces)
% 
% mb.applyForcesAndMoments (1)
% 
% pause (2);
% 
% mb.GetMotion ();
% 
% pos1 = mb.X(1)
% 
% mb.F (forces)
% 
% mb.applyForcesAndMoments (1)


clear mb;






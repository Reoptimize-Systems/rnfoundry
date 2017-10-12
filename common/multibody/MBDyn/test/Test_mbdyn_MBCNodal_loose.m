%% Test_mbdyn_MBCNodal


% start mbdyn 
delete ('output.*');

mbdpath = fullfile (getmfilepath ('Test_mbdyn_MBCNodal'), 'socket_loose.mbd');
% start mbdyn 

[status, cmdout] = system ( sprintf ('export LD_LIBRARY_PATH="" ; mbdyn -f "%s" -o "%s" > "%s" 2>&1 &', ...
                    mbdpath, ...
                    mbdpath(1:end-4), ...
                    [mbdpath(1:end-4), '.txt'] ...
                                     ) ...
                           );
                       
[status, cmdout] = system (sprintf ('mbdyn -f "%s" -o output > output.txt 2>&1 &', mbdpath))

% wait a few seconds for mbdyn to initialise
pause (3);

%% DataAndNext false

clear mb

mb = mbdyn.mint.MBCNodal ();

mb.Initialize ( 'local', '/tmp/mbdyn.sock', ...
                'NNodes', 2, ...
                'UseLabels', false, ...
                'Verbose', true, ...
                'DataAndNext', false );

nnodes = mb.GetNodes ()


maxiter = 3; % to match example program 

forces = [ 0., 0., 0.;
           0., 0., 1. ].';
       
status = mb.GetMotion ();

pos = mb.NodePositions ()

mb.F (forces)

mb.applyForcesAndMoments (false)
mb.applyForcesAndMoments (true)

ind = 2;
time = 1e-3;
status = 0;
while status == 0
    
    time = time + 1e-3;
    
    status = mb.GetMotion ();
    
    if status ~= 0
        pos = mb.NodePositions ()
        continue;
    end
    
    fprintf (1, 'final status: %d, ind: %d\n', status, ind);

    pos = mb.NodePositions ()
    
    rot = mb.GetRot()
    
%     for ind = 1:nnodes
%     
%         fprintf (1, 'Node %d has label %d\n', ind, mb.KinematicsLabel (ind));
%     
%     end

    if time > 0.25
        forces = [ 0., 0., 0.;
                   0., 1., 0. ].';
    end

    % set the forces
    mb.F (forces)

    mb.applyForcesAndMoments (false)
    
%     status = mb.GetMotion ();
    
%     pos = mb.NodePositions ()
    
%     mb.F (forces)

    mb.applyForcesAndMoments (true)
    
    ind = ind + 1;

end

clear mb;

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

%%


mbout = mbdyn.postproc ();

mbout.loadResultsFromFiles ( mbdpath );

mbout.plotNodeTrajectories ()

mbout.animate ()






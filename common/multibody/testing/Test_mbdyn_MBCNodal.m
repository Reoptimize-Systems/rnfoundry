%% Test_mbdyn_MBCNodal

% start mbdyn 
delete ('output.*');

if ispc
    mbdpath = fullfile (getmfilepath ('Test_mbdyn_MBCNodal'), 'socket_inet.mbd');
    precmd = '';
    mbdyncmd = ['"', fullfile(getmfilepath ('mexmbdyn_setup'), 'mbdyn_win64', 'bin', 'mbdyn.exe'), '"'];
else
    mbdpath = fullfile (getmfilepath ('Test_mbdyn_MBCNodal'), 'socket.mbd');
    precmd = 'export LD_LIBRARY_PATH="" ;';
    mbdyncmd = 'mbdyn';
end

[status, cmdout] = system ( sprintf ('%s %s -f "%s" -o "%s" > "%s" 2>&1 &', ...
                    precmd, ...
                    mbdyncmd, ...
                    mbdpath, ...
                    mbdpath(1:end-4), ...
                    [mbdpath(1:end-4), '.txt'] ...
                                     ) ...
                           );
                       
% [status, cmdout] = system (sprintf ('mbdyn -f "%s" -o output > output.txt 2>&1 &', mbdpath))

% wait a few seconds for mbdyn to initialise
% pause (3);

%%

if ispc
    mb = mbdyn.mint.MBCNodal ( 'CommMethod', 'inet socket', ...
                    'Host', '127.0.0.1', ...
                    'HostPort', 5500, ...
                    'NNodes', 2, ...
                    'UseLabels', false, ...
                    'MBDynInputFile', mbdpath);
else
    mb = mbdyn.mint.MBCNodal (  'CommMethod', 'local socket', ...
                    'Path', '/tmp/mbdyn.sock', ...
                    'NNodes', 2, ...
                    'UseLabels', false, ...
                    'MBDynInputFile', mbdpath);
end

mb.start ('Verbose', true);

nnodes = mb.GetNodes ();

maxiter = 3; % to match example program 

forces = [ 0., 0., 0.;
           0., 0., 0.1 ].';

status = 0;
while status == 0
    
    for iter = (maxiter-1):-1:1

%         disp(iter)
        status = mb.GetMotion ();
        fprintf (1, 'iter: %d, iterstatus: %d\n', iter, status);

%         pos1 = mb.X(1)
        
        % set the forces
        mb.F (forces)

        mb.applyForcesAndMoments (false)

    end
    
    status = mb.GetMotion ();
    fprintf (1, 'final status: %d\n', status);

    pos1 = mb.X(1)
    
    rot = mb.GetRot()
    
%     for ind = 1:nnodes
%     
%         fprintf (1, 'Node %d has label %d\n', ind, mb.KinematicsLabel (ind));
%     
%     end

    % set the forces
    mb.F (forces)

    mb.applyForcesAndMoments (true)

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


mbout = mbdyn.postproc (mbdpath);

% mbout.loadResultsFromFiles ( mbdpath );

mbout.plotNodeTrajectories ()








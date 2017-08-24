%% Test_mbdyn_MBCNodalSharedMem

% start mbdyn 
delete ('output.*');

devel = true;
    
if ispc
    mbdpath = fullfile (getmfilepath ('Test_mbdyn_MBCNodalSharedMem'), 'sharedmem.mbd');
    precmd = '';
    mbdyncmd = ['"', fullfile(getmfilepath ('mexmbdyn_setup'), 'mbdyn_win64', 'bin', 'mbdyn.exe'), '"'];
else
    mbdpath = fullfile (getmfilepath ('Test_mbdyn_MBCNodalSharedMem'), 'sharedmem.mbd');
    precmd = '';
    
    if devel
        mbdyncmd = '/opt/bin/mbdyn';
    else
        mbdyncmd = 'mbdyn';
    end
end

%%

% create the communicator object
mb = mbdyn.mint.MBCNodal ('CommMethod', 'shared memory', ...
                          'SharedMemoryName', 'shared_mem_location_name_43214321412', ...
                          'NNodes', 2, ...
                          'UseLabels', false, ...
                          'MBDynInputFile', mbdpath, ...
                          'MBDynExecutable', mbdyncmd);
                      
mb.start ('Verbosity', 3, 'MBDynStartWaitTime', 2);

nnodes = mb.GetNodes ()


maxiter = 3; % to match example program 

forces = [ 0., 0., 0.;
           0., 0., 0.1 ].';

status = 0;
while status == 0
    
    for iter = (maxiter-1):-1:1

%         disp(iter)
        status = mb.GetMotion ();
        fprintf (1, 'iter: %d, iterstatus: %d\n', iter, status);
        
        
        if status == -1
            break;
        end

%         pos1 = mb.X(1)
        
        % set the forces
        mb.F (forces)

        mb.applyForcesAndMoments (false)

    end
    
    if status == -1
        continue;
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

fprintf (1, 'Simulation complete');
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








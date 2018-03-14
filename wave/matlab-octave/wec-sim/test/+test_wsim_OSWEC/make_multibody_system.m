function mbsys = make_multibody_system (waves, simu, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options)

    default_problem_options.ResidualTol = 1e-6;
    default_problem_options.MaxIterations = 20;
    default_problem_options.Output = {'iterations', 'bailout'}; % , 'solution', 'jacobian matrix', 'matrix condition number', 'solver condition number'
    default_problem_options.NonLinearSolver = mbdyn.pre.newtonRaphsonSolver ();
    default_problem_options.LinearSolver = mbdyn.pre.linearSolver ('umfpack');
    default_problem_options.SteppingMethod = {};
        
    if nargin < 6
        problem_options = struct ();
    end
        
    problem_options = parseoptions (default_problem_options, problem_options);
    
    
    gref = mbdyn.pre.globalref;

    ref_seabed = mbdyn.pre.reference ( [0;0;-10], [], [], [], 'Parent', gref); 

    % Flap
    flap_node = hydro_mbnodes{1};

    flap_mb_body = hydro_mbbodies{1};

%     flap_mb_body.setSize (20, 20, 5);

    % make it a yellow colour
    flap_mb_body.setColour ([247, 244, 17]./255);

    % Base
    base_node = hydro_mbnodes{2};

    base_mb_body = hydro_mbbodies{2};

%     base_mb_body.setSize (6, 6, 38);

    % make it a grey colour
    base_mb_body.setColour ([181, 181, 181]./255);

    % clamp the base node to the global frame  so we can use to attach a
    % pin to the flap
    % clamp it
    jclamp = mbdyn.pre.clamp (base_node, 'node', 'node');

    % add the joints

    abs_hinge_pos = [0;0;-8.9];
    om = mbdyn.pre.orientmat ('2vectors', struct ('ia', 1, 'vecA', [1., 0., 0.], 'ib', 2, 'vecB', [0., 0., 1.]));
    jpin = mbdyn.pre.revolutePin ( flap_node, abs_hinge_pos, abs_hinge_pos, ...
                                   'PinOrientation', om);

    if ispc
        % add the socket forces
        communicator = mbdyn.pre.socketCommunicator ('Port', 5500, ...
                                                     'Coupling', 'tight', ...
                                                     'Create', 'yes', ...
                                                     'SleepTime', 0, ...
                                                     'SendAfterPredict', 'no');
    else
        % add the socket forces
        communicator = mbdyn.pre.socketCommunicator ('Path', '/tmp/mbdyn.sock', ...
                                                     'Coupling', 'tight', ...
                                                     'Create', 'yes', ...
                                                     'SleepTime', 0, ...
                                                     'SendAfterPredict', 'no');
    end

    socket_force = mbdyn.pre.externalStructuralForce ({flap_node, base_node}, [], communicator, ...
                                                      'Labels', 'no', ...
                                                      'Sorted', 'yes', ...
                                                      'Orientation', 'orientation matrix' ... , ...
                                                      ... 'Echo', fullfile (simu.caseDir, 'output', 'socket_echo.txt') ...
                                                      );

    prob = mbdyn.pre.initialValueProblem (simu.startTime, simu.endTime, simu.dt, ...
                                    'ResidualTol', problem_options.ResidualTol, ...
                                    'MaxIterations', problem_options.MaxIterations, ...
                                    'Output', problem_options.Output, ...
                                    'NonlinearSolver', problem_options.NonLinearSolver, ...
                                    'LinearSolver', problem_options.LinearSolver, ...
                                    'Method', problem_options.SteppingMethod ); 

    % assemble the system
    mbsys = mbdyn.pre.system ( prob, ...
                               'Nodes', hydro_mbnodes, ...
                               'Elements', [{flap_mb_body, base_mb_body, jclamp, jpin, socket_force}, hydro_mbelements], ...
                               'DefaultOrientation', 'orientation matrix', ...
                               'DefaultOutput', {'none'} );

end
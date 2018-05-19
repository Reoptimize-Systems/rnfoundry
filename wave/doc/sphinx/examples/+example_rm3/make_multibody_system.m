function [mbsys, initptodpos] = make_multibody_system (waves, simu, hydro_mbnodes, hydro_mbbodies, hydro_mbelements, problem_options)

    default_problem_options.ResidualTol = 1e-6;
    default_problem_options.MaxIterations = 20;
    default_problem_options.Output = {'iterations', 'bailout'}; % , 'solution', 'jacobian matrix', 'matrix condition number', 'solver condition number'
    default_problem_options.NonLinearSolver = mbdyn.pre.newtonRaphsonSolver ();
    default_problem_options.LinearSolver = mbdyn.pre.linearSolver ('umfpack');
    default_problem_options.SteppingMethod = {};
    default_problem_options.DefaultElementOutput = {'all'};
        
    if nargin < 6
        problem_options = struct ();
    end
        
    problem_options = parseoptions (default_problem_options, problem_options);
    
    % make reference to the gobal frame, other references we create will be
    % relative to this one. The global origin is the mean water surface
    % level in this case
    gref = mbdyn.pre.globalref;

    % % Get the float node (created previously by the hydro system and
    % supplied as an input to this function).
    float_node = hydro_mbnodes{1};

    % similarly get the float body
    float_mb_body = hydro_mbbodies{1};

    % make it a yellow colour for visualisation purposes
    float_mb_body.setColour ([247, 244, 17]./255);

    % Get the spar node (created previously by the hydro system and
    % supplied as an input to this function).
    spar_node = hydro_mbnodes{2};

    % similarly get the spar body 
    spar_mb_body = hydro_mbbodies{2};

    % make it a grey colour
    spar_mb_body.setColour ([181, 181, 181]./255);

    % next create a reference for the sea bed we will put a node on the sea
    % bed which we will use to add some constraints on the other nodes,
    % actually we could upt this node anywhere, but the sea bed is as good
    % a place as any
    ref_seabed = mbdyn.pre.reference ( [0;0;-waves.waterDepth], [], [], [], 'Parent', gref); 
    
    % create a node to clamp to the global frame which we can use to attach
    % a planar joint (actually made up of two joints) to the spar. The
    % planar joint will force the spar to move in a plane attached to the
    % clamped node. This restricts the motion of the whole system to the
    % X-Z plane (this is to match the RM3 simulation from the original
    % WECSim and is not essential for any system).
    %
    % Since  we are going to clamp the node, we can make it static (getting rid
    % of the degrees of freedom associated with it)
    clamped_node = mbdyn.pre.structuralNode6dof ('static', 'AbsolutePosition', ref_seabed.pos);
    
    % clamp it with a clamp joint
    jclamp = mbdyn.pre.clamp (clamped_node, 'node', 'node');

    % apply the inPlane constraint so spar can only move in X-Z axis
    j3 = mbdyn.pre.inPlane (  clamped_node, spar_node, [0;1;0], ...
                              'RelativeOffset', 'null', ...
                              'RelativeOffsetReference', 'other node');
                          
    % add a prismatic joint (keeps the orientation of two nodes fixed, i.e.
    % they won't rotate relative to each other)
    j1 = mbdyn.pre.prismatic (float_node, spar_node);

    
    % The next section creates the pto mechanism constraints
    
    % create an orientation 
    om = mbdyn.pre.orientmat ( '2vectors', ...
                               struct ( 'ia', 1, ...
                                        'vecA', [1., 0., 0.], ...
                                        'ib', 3, ...
                                        'vecB', [0., 0., 1.] ) ...
                             );
                           
    j2 = mbdyn.pre.inline ( float_node, spar_node, ...
                            'RelativeLinePosition', [0;0;0], ...
                            'RelativeOrientation', om, ...
                            'OrientationReference', 'global' );


	% create an orientation with axis 3 pointing along the global axis 2
	% (the gloabl Y axis)
    Yax_orient = mbdyn.pre.orientmat ( '2vectors', ...
                                        struct ( 'ia', 1, ...
                                                 'vecA', [1,0,0], ...
                                                 'ib', 2, ...
                                                 'vecB', [0,0,1] ) ...
                                      );
    
	% use it in a revolute rotation constraint to only allow the WEC to
	% rotate in the X-Z plane. The revolute rotation element allows roation
	% only about axis 3 of the supplied orientaion
    j4 = mbdyn.pre.revoluteRotation ( spar_node, clamped_node, ...
                                     'RelativeOffset1', 'null', ...
                                     'RelativeOffset1Reference', 'node', ...
                                     'RelativeOffset2', 'null', ...
                                     'RelativeOffset2Reference', 'other node', ...
                                     'RelativeOrientation1', Yax_orient, ...
                                     'RelativeOrientation2', Yax_orient );


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

    socket_force = mbdyn.pre.externalStructuralForce ({float_node, spar_node}, [], communicator, ...
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
                               'Nodes', [hydro_mbnodes, {clamped_node}], ...
                               'Elements', [{float_mb_body, spar_mb_body, jclamp, j1, j2, j3, j4 socket_force}, hydro_mbelements], ...
                               'DefaultOrientation', 'orientation matrix', ...
                               'DefaultOutput', problem_options.DefaultElementOutput );
                           
	initptodpos = float_node.absolutePosition - spar_node.absolutePosition; 

end
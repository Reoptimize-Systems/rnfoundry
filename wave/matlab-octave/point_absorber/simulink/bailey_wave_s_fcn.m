function bailey_wave_s_fcn (block)
% A level-2 s-function representation of a three phase generator emf.  This
% has been setup using a standard level-2 s-function template.
%
% See also: rotgensfcn_common.m

% Created Richard Crozier 2013

    % Call setup function
    setup(block);

end

function setup(block)
% Sets up the s-function block's basic characteristics such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
    
    % Register number of ports
    % Define the number of input ports (heave and surge force on point absorber)
    block.NumInputPorts  = 2;

    % Define the number of output ports (position in heave and surge and velocity in heave and surge, then other data multiplexed)
    block.NumOutputPorts = 5;
    
    % Setup port properties to be inherited or dynamic
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    % Define Block Inputs
    % First input is heave force
    block.InputPort(1).Dimensions  = 1;
    block.InputPort(1).DatatypeID  = 0;  % double
    block.InputPort(1).Complexity  = 'Real';
    block.InputPort(1).DirectFeedthrough = true;

    % Second input is surge force
    block.InputPort(2).Dimensions  = 1;
    block.InputPort(2).DatatypeID  = 0;  % double
    block.InputPort(2).Complexity  = 'Real';
    block.InputPort(2).DirectFeedthrough = true;
    
    % Define Block Outputs
    % First output is the position in heave
    block.OutputPort(1).Dimensions  = 1;
    block.OutputPort(1).DatatypeID  = 0; % double
    block.OutputPort(1).Complexity  = 'Real';

    % Second output is the position in surge
    block.OutputPort(2).Dimensions  = 1;
    block.OutputPort(2).DatatypeID  = 0; % double
    block.OutputPort(2).Complexity  = 'Real';
       
    % Third output is the velocity in heave
    block.OutputPort(3).Dimensions  = 1;
    block.OutputPort(3).DatatypeID  = 0; % double
    block.OutputPort(3).Complexity  = 'Real';
    
    % Fourth output is the velocity in surge
    block.OutputPort(4).Dimensions  = 1;
    block.OutputPort(4).DatatypeID  = 0; % double
    block.OutputPort(4).Complexity  = 'Real';
    
    % Fifth output is various quantities multiplexed together
    block.OutputPort(5).Dimensions  = 8;
    block.OutputPort(5).DatatypeID  = 0; % double
    block.OutputPort(5).Complexity  = 'Real';
    
    
    % Register Dialog parameters
    block.NumDialogPrms     = 3;
%     block.NumDialogPrms     = 1;
%     block.DialogPrmsTunable = {'Nontunable'};

    % Register sample times ([0 offset] = Continuous sample time)
    block.SampleTimes = [0 0];
    
    % Continuous states (derivatives), the position and velocity in heave
    % and surge and the radiation states
%     block.NumContStates = 24;
    NRadiationCoefs = block.DialogPrm(3).Data;
    block.NumContStates = 4 + 2*NRadiationCoefs;

    % Register all relevant block methods
%     block.RegBlockMethod('CheckParameters',      @CheckPrms);
    block.RegBlockMethod('PostPropagationSetup',@DoPostPropSetup);
    block.RegBlockMethod('Start', @Start);
    block.RegBlockMethod('Outputs', @Outputs);
    block.RegBlockMethod('InitializeConditions', @InitConditions);
    block.RegBlockMethod('Derivatives', @Derivatives);
    block.RegBlockMethod('Terminate', @Terminate);
    block.RegBlockMethod('SetInputPortSamplingMode',@SetInputPortSamplingMode);
    block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);
    
end

function DoPostPropSetup(block)

    % set the number of radiation coefficients
    NRadiationCoefs = block.DialogPrm(3).Data;
    block.NumContStates = 4 + 2*NRadiationCoefs;
  
end

function Start(block)
% check and initialises block just prior to simulation

    % get the buoy choice
    buoy = block.DialogPrm(1).Data;
    
    % get the sea pamaters structure (should have been generated with
    % seasetup.m)
    buoysimoptions.SeaParameters = block.DialogPrm(2).Data;
    
    % get the desired number of radiation coefficients
    buoysimoptions.NRadiationCoefs = block.DialogPrm(3).Data;

    % set up the buoy in readiness for simulation
    buoysimoptions = buoysimsetup (buoy, buoysimoptions);

    % put the data into UserData for the block
    set(block.BlockHandle, 'UserData', struct ('buoysimoptions', buoysimoptions));

end

function SetInputPortSamplingMode(block, idx, fd)
% Set the port sampling modes
    
    % Input port
    block.InputPort(idx).SamplingMode = fd;

    % Output ports
    for i = 1:block.NumOutputPorts
        block.OutputPort(i).SamplingMode = fd;
    end
    
end

function SetInpPortDims(block, idx, di)
% Set input port dimensions
    block.InputPort(idx).Dimensions = di;
end

function Outputs(block)
% Outputs - called to generate block outputs in simulation step

    % get the data structure containing design and simparams from the block
    % UserData parameter
    data = get(block.BlockHandle, 'UserData');

    [ bouyancy_force, ...
      excitation_force_heave, ...
      excitation_force_surge, ...
      radiation_force_heave, ...
      radiation_force_surge, ...
      FBDh, ...
      FBDs, ...
      wave_height ] = buoyodeforces ( block.CurrentTime, ...
                                      block.ContStates.Data(5:end), ...
                                      block.ContStates.Data(1), ...
                                      block.ContStates.Data(2), ...
                                      block.ContStates.Data(4), ...
                                      data.buoysimoptions );

    % put the output data in the block output port data fields 
    
    % heave position
    block.OutputPort(1).Data = block.ContStates.Data(1);
    block.OutputPort(2).Data = block.ContStates.Data(2);
    block.OutputPort(3).Data = block.ContStates.Data(3);
    block.OutputPort(4).Data = block.ContStates.Data(4);

    % other outputs
    block.OutputPort(5).Data = [ bouyancy_force, ...
                                 excitation_force_heave, ...
                                 excitation_force_surge, ...
                                 radiation_force_heave, ...
                                 radiation_force_surge, ...
                                 FBDh, ...
                                 FBDs, ...
                                 wave_height ];
                                 
    
end


function InitConditions(block)
% construct the initial conditions of the system
    
    data = get(block.BlockHandle, 'UserData');
    
    block.ContStates.Data = [ ...
        data.buoysimoptions.ODESim.SolutionComponents.BuoyPositionHeave.InitialConditions, ...
        data.buoysimoptions.ODESim.SolutionComponents.BuoyVelocityHeave.InitialConditions, ...
        data.buoysimoptions.ODESim.SolutionComponents.BuoyPositionSurge.InitialConditions, ...
        data.buoysimoptions.ODESim.SolutionComponents.BuoyVelocitySurge.InitialConditions, ...
        data.buoysimoptions.ODESim.SolutionComponents.BuoyRadiationHeave.InitialConditions, ...
        data.buoysimoptions.ODESim.SolutionComponents.BuoyRadiationSurge.InitialConditions ...
                            ];
    
end

function Derivatives(block)
% Derivatives - called to generate block derivatives at each simulation
% step

    data = get(block.BlockHandle, 'UserData');
    
    % calculate the derivatives, first four are velocity and acceleration
    % in heave and surge
    dx = buoyodesim ( block.CurrentTime, ...
                      block.ContStates.Data, ...
                      data.buoysimoptions, ...
                      [block.InputPort(1).Data, block.InputPort(2).Data] );
    
    % get the derivatives
    block.Derivatives.Data = dx;

end
    
function Terminate(block) %#ok<INUSD>
% Terminate block at the end of simulation
end
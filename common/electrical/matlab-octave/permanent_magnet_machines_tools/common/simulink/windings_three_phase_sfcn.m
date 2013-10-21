function rotgen_three_phase_sfcn(block)
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
    % Define the number of input ports (position, omega, and currents)
    block.NumInputPorts  = 7;

    % Define the number of output ports (torque and phase EMFs)
    block.NumOutputPorts = 1;
    
    % Setup port properties to be inherited or dynamic
%     block.SetPreCompInpPortInfoToDynamic;
%     block.SetPreCompOutPortInfoToDynamic;

    % Define Block Inputs;

    % phase A resistance
    block.InputPort(1).Dimensions  = 1;
    block.InputPort(1).DatatypeID  = 0;  % double
    block.InputPort(1).Complexity  = 'Real';
%     block.InputPort(1).DirectFeedthrough = true;
    
    % phase B resistance
    block.InputPort(2).Dimensions  = 1;
    block.InputPort(2).DatatypeID  = 0;  % double
    block.InputPort(2).Complexity  = 'Real';
%     block.InputPort(2).DirectFeedthrough = true;
    
    % phase C resistance
    block.InputPort(3).Dimensions  = 1;
    block.InputPort(3).DatatypeID  = 0;  % double
    block.InputPort(3).Complexity  = 'Real';
%     block.InputPort(3).DirectFeedthrough = true;
    
    % phase A inductance
    block.InputPort(4).Dimensions  = 1;
    block.InputPort(4).DatatypeID  = 0;  % double
    block.InputPort(4).Complexity  = 'Real';
%     block.InputPort(4).DirectFeedthrough = true;
    
    % phase B inductance
    block.InputPort(5).Dimensions  = 1;
    block.InputPort(5).DatatypeID  = 0;  % double
    block.InputPort(5).Complexity  = 'Real';
%     block.InputPort(5).DirectFeedthrough = true;
    
    % phase C inductance
    block.InputPort(6).Dimensions  = 1;
    block.InputPort(6).DatatypeID  = 0;  % double
    block.InputPort(6).Complexity  = 'Real';
%     block.InputPort(6).DirectFeedthrough = true;
    
    % mutual inductance
    block.InputPort(7).Dimensions  = 1;
    block.InputPort(7).DatatypeID  = 0;  % double
    block.InputPort(7).Complexity  = 'Real';
%     block.InputPort(7).DirectFeedthrough = true;
    
    % output is the calculated mutual resistance
    block.OutputPort(1).Dimensions  = 1;
    block.OutputPort(1).DatatypeID  = 0; % double
    block.OutputPort(1).Complexity  = 'Real';
    
    % Register Dialog parameters
    block.NumDialogPrms     = 0;
%     block.NumDialogPrms     = 1;
%     block.DialogPrmsTunable = {'Nontunable'};

    % Register sample times ([0 offset] = Continuous sample time)
    block.SampleTimes = [0 0];

    % Register all relevant block methods
%     block.RegBlockMethod('CheckParameters',      @CheckPrms);
    block.RegBlockMethod('Outputs', @Outputs);
    block.RegBlockMethod('Terminate', @Terminate);
    block.RegBlockMethod('SetInputPortSamplingMode',@SetInputPortSamplingMode);
    block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);
    
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

    % get the winding parameters data from the input ports
    phaseAmainresistance = block.InputPort(1).Data;
    phaseBmainresistance = block.InputPort(2).Data;
    phaseCmainresistance = block.InputPort(3).Data;
    
    phaseAmaininductance = block.InputPort(4).Data;
    phaseBmaininductance = block.InputPort(5).Data;
    phaseCmaininductance = block.InputPort(6).Data;
    
    mutualinductance   = block.InputPort(7).Data;
    mutualresistance = 1e6 * mean([phaseAmainresistance, phaseBmainresistance, phaseCmainresistance]);
    
    block.OutputPort(1).Data = mutualresistance;
    
    % apply the windings data to the three-phase mutual windings block
    set_param('Windings', 'SelfImpedance1', [phaseAmainresistance, phaseAmaininductance]);
    set_param('Windings', 'SelfImpedance2', [phaseBmainresistance, phaseBmaininductance]);
    set_param('Windings', 'SelfImpedance3', [phaseCmainresistance, phaseCmaininductance]);
    set_param('Windings', 'MutualImpedance', [mutualresistance, mutualinductance]);
    
end
    
function Terminate(block) %#ok<INUSD>
% Terminate block at the end of simulation
end
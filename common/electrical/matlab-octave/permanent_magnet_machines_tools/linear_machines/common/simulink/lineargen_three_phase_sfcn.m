function lineargen_three_phase_sfcn(block)
% A level-2 s-function representation of a three phase generator emf.  This
% has been setup using a standard level-2 s-function template.
%
% See also: lineargensfcn_common.m

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
    block.NumInputPorts  = 5;

    % Define the number of output ports (torque and phase EMFs)
    block.NumOutputPorts = 5;
    
    % Setup port properties to be inherited or dynamic
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;

    % Define Block Inputs
    % First input is position
    block.InputPort(1).Dimensions  = 1;
    block.InputPort(1).DatatypeID  = 0;  % double
    block.InputPort(1).Complexity  = 'Real';
    block.InputPort(1).DirectFeedthrough = true;

    % Second input is velocity
    block.InputPort(2).Dimensions  = 1;
    block.InputPort(2).DatatypeID  = 0;  % double
    block.InputPort(2).Complexity  = 'Real';
    block.InputPort(2).DirectFeedthrough = true;

    % Third input is the currents in phase A
    block.InputPort(3).Dimensions  = 1;
    block.InputPort(3).DatatypeID  = 0;  % double
    block.InputPort(3).Complexity  = 'Real';
    block.InputPort(3).DirectFeedthrough = false;
    
    % Fourth input is the currents in phase B
    block.InputPort(4).Dimensions  = 1;
    block.InputPort(4).DatatypeID  = 0;  % double
    block.InputPort(4).Complexity  = 'Real';
    block.InputPort(4).DirectFeedthrough = false;
    
    % Fifth input is the currents in phase C
    block.InputPort(5).Dimensions  = 1;
    block.InputPort(5).DatatypeID  = 0;  % double
    block.InputPort(5).Complexity  = 'Real';
    block.InputPort(5).DirectFeedthrough = false;
    
    % Define Block Outputs
    % First output is the electromagnetic torque provided by the generator
    block.OutputPort(1).Dimensions  = 1;
    block.OutputPort(1).DatatypeID  = 0; % double
    block.OutputPort(1).Complexity  = 'Real';

    % Second output is the EMFs in phase A
    block.OutputPort(2).Dimensions  = 1;
    block.OutputPort(2).DatatypeID  = 0; % double
    block.OutputPort(2).Complexity  = 'Real';
       
    % Third output is the EMFs in phase B
    block.OutputPort(3).Dimensions  = 1;
    block.OutputPort(3).DatatypeID  = 0; % double
    block.OutputPort(3).Complexity  = 'Real';
    
    % Fourth output is the EMFs in phase C
    block.OutputPort(4).Dimensions  = 1;
    block.OutputPort(4).DatatypeID  = 0; % double
    block.OutputPort(4).Complexity  = 'Real';
    
    % Fifth output is various quantities multiplexed together
    block.OutputPort(5).Dimensions  = 12;
    block.OutputPort(5).DatatypeID  = 0; % double
    block.OutputPort(5).Complexity  = 'Real';
    
    
    % Register Dialog parameters
    block.NumDialogPrms     = 0;
%     block.NumDialogPrms     = 1;
%     block.DialogPrmsTunable = {'Nontunable'};

    % Register sample times ([0 offset] = Continuous sample time)
    block.SampleTimes = [0 0];
    
    % Continuous states (derivatives), the angular accleration
%     block.NumContStates = 1;
    block.NumContStates = 0;

    % Register all relevant block methods
%     block.RegBlockMethod('CheckParameters',      @CheckPrms);
    block.RegBlockMethod('Start', @Start);
    block.RegBlockMethod('Outputs', @Outputs);
%     block.RegBlockMethod('Derivatives', @Derivatives);
    block.RegBlockMethod('Terminate', @Terminate);
    block.RegBlockMethod('SetInputPortSamplingMode',@SetInputPortSamplingMode);
    block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);
    
end

function Start(block)
% check and initialises block just prior to simulation

    % get the UserData for the block
    data = get (block.BlockHandle, 'UserData');
    
    if isempty(data)
        % get the data from the parameters

        % get the value of the drop down to see if we are getting design data
        % from file or workspace
        file_or_workspace = block.DialogPrm(1).Data;

        if strcmp (file_or_workspace, 'From File')

            % get the file name
            design_data_file = block.DialogPrm(2).Data;

            % check it exists
            if exist (design_data_file, 'file')
                % check if right data is in file
                try
                    load ('design_data_file', 'design', 'simoptions');
                catch ME
                    error ('RENEWNET:lineargen_three_phase_sfcn:badfile', ...
                    'Loading design data from file\n: %s\nfailed with the following error message:\n"%s"\n', ...
                        design_data_file, ME.message);
                end

            else
                error ('RENEWNET:lineargen_three_phase_sfcn:badfile', ...
                    'Design data file\n: %s\ndoes not appear to exist', design_data_file);
            end

        elseif strcmp (file_or_workspace, 'From Workspace')
            % get the variable names 
            design = block.DialogPrm(3).Data;
            simoptions = block.DialogPrm(4).Data;
        else
            error ('RENEWNET:lineargen_three_phase_sfcn:baddatatype', ...
                'Invalid Data Type popup value')
        end

        % do some checks on the data

        if design.Phases ~= 3
            error('Design must have three Phases.')
        end

        if ~isfield(design, 'slm_fluxlinkage')
            error('Design data does not appear to be complete (no flux linage info)')
        end

        % put the data into UserData for the block
        set(block.BlockHandle, 'UserData', struct ('design', design, ...
                                                   'simoptions', simoptions))
                                           
    end

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

    % get the simulation data from the input ports
    xR                  = block.InputPort(1).Data;
    vR                  = block.InputPort(2).Data;
    phasecurrent(1,1)   = block.InputPort(3).Data;
    phasecurrent(2,1)   = block.InputPort(4).Data;
    phasecurrent(3,1)   = block.InputPort(5).Data;

    [EMF, force, flux, data.design] = rotgensfcn_common(data.design, ...
                                      data.simoptions,...
                                      xR,...
                                      vR, ...
                                      phasecurrent);

    % put the output data in the block output port data fields 
    
    % torque
    block.OutputPort(1).Data = force;
    % emfs
    block.OutputPort(2).Data = EMF(1);
    block.OutputPort(3).Data = EMF(2);
    block.OutputPort(4).Data = EMF(3);
    % phase resistances etc
    block.OutputPort(5).Data = [ data.design.RPhase(1,1), ... % the phase resistance of phase A
                                 data.design.RPhase(2,2), ... % the phase resistance of phase B
                                 data.design.RPhase(3,3), ... % the phase resistance of phase C
                                 flux', ...  % the flux linkage in the Phases
                                 vR / (2 * data.design.PoleWidth), ... % the electrical frequency
                                 vR, ...
                                 2 * pi * xR / (2 * data.design.PoleWidth), ... % the electrical angle
                                 force, ...
                                 data.design.RotorMomentOfInertia, ...
                                 data.design.FluxLinkageRms ];
    
end

% function Derivatives(block)
% % Derivatives - called to generate block derivatives at each simulation
% % step
% 
%     data = get(block.BlockHandle, 'UserData');
%     
%     % the angular acceleration, the derivative of 
%     block.Derivatives.Data = 
%     block.InputPort(1).Data;
%     block.ContStates.Data
% 
% end
    
function Terminate(block) %#ok<INUSD>
% Terminate block at the end of simulation
end